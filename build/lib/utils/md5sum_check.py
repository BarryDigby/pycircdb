import csv
import hashlib
import os
from typing import Dict, List, Optional


def _file_md5sum(file_path: str) -> str:
    """Calculate the MD5 checksum for a single file."""
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def _load_expected_sums(md5sum_file: str) -> Dict[str, str]:
    """Load checksum manifest as {filename: md5sum}."""
    with open(md5sum_file, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        normalized = {name.lower().strip(): name for name in fieldnames}

        file_col = normalized.get("file") or normalized.get("filename")
        sum_col = normalized.get("md5sum") or normalized.get("md5")
        if not file_col or not sum_col:
            raise ValueError(
                "md5sum_file must include file/filename and md5sum/md5 columns"
            )

        expected = {}
        for row in reader:
            filename = (row.get(file_col) or "").strip()
            checksum = (row.get(sum_col) or "").strip().lower()
            if filename and checksum:
                expected[filename] = checksum

    if not expected:
        raise ValueError("md5sum_file contains no checksum entries")

    return expected


def _candidate_dirs(tmp_dir: str, dir_prefix: str) -> List[str]:
    """Return candidate tmp subdirectories sorted newest-first."""
    dirs = []
    for entry in os.listdir(tmp_dir):
        dir_path = os.path.join(tmp_dir, entry)
        if os.path.isdir(dir_path) and entry.startswith(dir_prefix):
            dirs.append(dir_path)
    dirs.sort(key=os.path.getmtime, reverse=True)
    return dirs


def _validate_dir(
    dir_path: str,
    expected: Dict[str, str],
    files_to_validate: List[str],
) -> List[str]:
    """Return validated file paths for a directory, or [] if any file fails."""
    valid_paths: List[str] = []

    for filename in files_to_validate:
        expected_sum = expected.get(filename)
        if not expected_sum:
            return []

        file_path = os.path.join(dir_path, filename)
        if not os.path.isfile(file_path):
            return []

        if _file_md5sum(file_path).lower() != expected_sum:
            return []

        valid_paths.append(file_path)

    return valid_paths


def check_sums(
    tmp_dir: str,
    md5sum_file: str,
    dir_prefix: str = "lookup_tables_",
    required_files: Optional[List[str]] = None,
) -> List[str]:
    """
    Validate parquet files in tmp subdirectories against a checksum manifest.

    Returns a list of valid parquet file paths when a fully valid directory is found.
    Returns an empty list when no valid directory exists.
    """
    expected = _load_expected_sums(md5sum_file)
    files_to_validate = required_files or list(expected.keys())

    if not os.path.isdir(tmp_dir):
        return []

    candidate_dirs = _candidate_dirs(tmp_dir=tmp_dir, dir_prefix=dir_prefix)

    for dir_path in candidate_dirs:
        valid_paths = _validate_dir(
            dir_path=dir_path,
            expected=expected,
            files_to_validate=files_to_validate,
        )
        if valid_paths:
            print(f"Using cached files from {dir_path}; all MD5 checks passed.")
            return valid_paths

    print("No cached lookup directory passed MD5 validation; fresh download required.")
    return []