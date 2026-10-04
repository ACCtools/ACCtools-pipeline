"""Bind newly generated assembly alignments to their actual source sequences.

This module does not attest existing PAFs retrospectively. Its separate index
cache binds reference content at index generation, including the mapper binary.
"""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


INDEX_SCHEMA = "SKYPE.reference_index_source.v1"
SOURCE_SCHEMA = "SKYPE.assembly_alignment_source.v1"


def content_signature(path):
    """Hash a stable regular input and retain its resolved path."""
    resolved = Path(path).resolve(strict=True)
    with resolved.open("rb") as handle:
        before = os.fstat(handle.fileno())
        digest = hashlib.sha256()
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
        after = os.fstat(handle.fileno())
    current = resolved.stat()
    identity = lambda stat: (stat.st_dev, stat.st_ino, stat.st_size,
                             stat.st_mtime_ns, stat.st_ctime_ns)
    if identity(before) != identity(after) or identity(after) != identity(current):
        raise RuntimeError(f"Input changed while computing its content binding: {path}")
    return {"path": str(resolved), "sha256": digest.hexdigest()}


def write_json_atomic(path, data):
    path = Path(path)
    fd, temporary = tempfile.mkstemp(prefix=".source-binding-", suffix=".json.tmp",
                                     dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(data, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def mapper_identity():
    binary = shutil.which("minimap2")
    if binary is None:
        raise FileNotFoundError("minimap2 is required to create assembly alignments")
    identity = content_signature(binary)
    version = subprocess.check_output([identity["path"], "--version"], text=True).strip()
    if identity != content_signature(binary):
        raise RuntimeError("minimap2 changed while obtaining its identity")
    return dict(identity, version=version)


def _valid_index_record(record, reference, mapper, index_path, preset):
    try:
        command = record["generation_command"]
        return (record["schema"] == INDEX_SCHEMA and record["status"] == "complete"
                and record["reference"] == reference
                and record["minimap2"] == mapper and record["preset"] == preset
                and record["index"] == content_signature(index_path)
                and command[0] == mapper["path"] and command[-1] == reference["path"]
                and command[1:3] == ["-x", preset] and "-d" in command
                and command[command.index("-d")+1] == record["generation_output"]["path"]
                and record["generation_output"]["sha256"] == record["index"]["sha256"])
    except (KeyError, IndexError, TypeError, ValueError, OSError):
        return False


def ensure_bound_reference_index(reference, cache_dir, thread, mapper=None):
    """Return an asm20 index whose generation was bound to reference content.

    Existing size/mtime-only indexes are left intact. They cannot establish this
    relationship, so the first bound run builds a distinct shared index.
    """
    preset = "asm20"
    source = content_signature(reference)
    mapper = mapper or mapper_identity()
    key_data = dict(reference=source, minimap2=mapper, preset=preset)
    key = hashlib.sha256(json.dumps(key_data, sort_keys=True).encode()).hexdigest()
    cache_dir = Path(cache_dir).resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)
    index = cache_dir / f"{Path(source['path']).name}.source-bound.{key[:24]}.mmi"
    manifest = Path(str(index) + ".source_binding.json")
    with Path(str(index) + ".lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if manifest.is_file() and index.is_file() and index.stat().st_size:
            try:
                record = json.loads(manifest.read_text())
            except (ValueError, OSError):
                record = {}
            if _valid_index_record(record, source, mapper, index, preset):
                return record
        # A failed build must not retain an earlier successful binding.
        manifest.unlink(missing_ok=True)
        with tempfile.TemporaryDirectory(prefix="source-bound-index-", dir=cache_dir) as temporary:
            output = Path(temporary) / "reference.mmi"
            command = [mapper["path"], "-x", preset, "-t", str(thread),
                       "-d", str(output), source["path"]]
            with Path(str(index) + ".log").open("w") as log:
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
            if not output.is_file() or not output.stat().st_size:
                raise RuntimeError("minimap2 completed without a nonempty reference index")
            if source != content_signature(reference):
                raise RuntimeError("Reference changed during index generation")
            if {k: mapper[k] for k in ("path", "sha256")} != content_signature(mapper["path"]):
                raise RuntimeError("minimap2 changed during index generation")
            generated_output = content_signature(output)
            os.replace(output, index)
            published_index = content_signature(index)
            if generated_output["sha256"] != published_index["sha256"]:
                raise RuntimeError("Index changed while publishing its completed output")
            record = dict(schema=INDEX_SCHEMA, status="complete", reference=source,
                          index=published_index, preset=preset, minimap2=mapper,
                          generation_command=command,
                          generation_output=generated_output,
                          output_publication="Successful temporary -d output atomically renamed to index.path")
            write_json_atomic(manifest, record)
        return record


def source_inputs(fasta, reference):
    return dict(fasta=content_signature(fasta), reference=content_signature(reference))


def publish_source_binding(path, before, fasta, reference, primary, alternate,
                           commands, index_binding, gap_extractor, generated_outputs):
    """Publish only after successful primary and alternate generation."""
    after = source_inputs(fasta, reference)
    if before != after:
        raise RuntimeError("Assembly or reference changed during alignment generation")
    mapper = index_binding["minimap2"]
    index_path = index_binding["index"]["path"]
    if not _valid_index_record(index_binding, before["reference"], mapper,
                               index_path, "asm20"):
        raise RuntimeError("Reference index binding changed during alignment generation")
    if {k: mapper[k] for k in ("path", "sha256")} != content_signature(mapper["path"]):
        raise RuntimeError("minimap2 changed during alignment generation")
    if gap_extractor != content_signature(gap_extractor["path"]):
        raise RuntimeError("Gap extraction code changed during alignment generation")
    mapping_commands = [command for command in commands if command[0] == mapper["path"]]
    if not mapping_commands or any(index_path not in command for command in mapping_commands):
        raise RuntimeError("Alignment commands do not use the bound reference index")
    current_outputs = dict(primary_paf=content_signature(primary),
                           alternate_paf=content_signature(alternate))
    if current_outputs != generated_outputs:
        raise RuntimeError("Generated primary or alternate PAF changed before provenance publication")
    record = dict(schema=SOURCE_SCHEMA, status="complete",
                  producer_stage="raw_and_alternate_alignment_generation",
                  inputs_before=before, inputs_after=after,
                  outputs=current_outputs, outputs_at_generation=generated_outputs,
                  generation_commands=commands, reference_index_binding=index_binding,
                  gap_extractor=gap_extractor,
                  limitations=["Content provenance is not evidence of alignment accuracy.",
                               "An empty primary alignment deterministically produces an empty alternate PAF."])
    write_json_atomic(path, record)
    return record
