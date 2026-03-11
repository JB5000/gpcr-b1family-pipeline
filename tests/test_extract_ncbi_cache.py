#!/usr/bin/env python3
"""
Unit tests for extract_ncbi_cache.py
Tests cover the pure (offline) helper functions:
  - sanitize_name
  - find_input_files (with temp directories)
  - parse_phmmer_file
  - print_progress
  - resolve_cache_out
"""

import os
import tempfile
import time
import unittest

from extract_ncbi_cache import (
    find_input_files,
    parse_phmmer_file,
    print_progress,
    resolve_cache_out,
    sanitize_name,
)


class TestSanitizeName(unittest.TestCase):
    def test_alphanumeric_unchanged(self):
        self.assertEqual(sanitize_name("CALCR"), "CALCR")

    def test_spaces_replaced(self):
        result = sanitize_name("my input file")
        self.assertNotIn(" ", result)

    def test_special_chars_replaced(self):
        result = sanitize_name("file-name.txt!")
        self.assertTrue(result.isidentifier() or result.replace("_", "").isalnum())

    def test_empty_fallback(self):
        self.assertEqual(sanitize_name(""), "output")

    def test_only_symbols_fallback(self):
        self.assertEqual(sanitize_name("---"), "output")


class TestFindInputFiles(unittest.TestCase):
    def test_finds_txt_files(self):
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, "sample.txt"), "w").close()
            open(os.path.join(d, "sample2.tsv"), "w").close()
            open(os.path.join(d, "ignore.csv"), "w").close()
            files = find_input_files(d)
            basenames = [os.path.basename(f) for f in files]
            self.assertIn("sample.txt", basenames)
            self.assertIn("sample2.tsv", basenames)
            self.assertNotIn("ignore.csv", basenames)

    def test_raises_on_missing_dir(self):
        with self.assertRaises(SystemExit):
            find_input_files("/nonexistent_dir_xyz_123")

    def test_raises_on_empty_dir(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(SystemExit):
                find_input_files(d)


class TestParsePhmmerFile(unittest.TestCase):
    def _write_temp(self, content: str) -> str:
        fd, path = tempfile.mkstemp(suffix=".txt")
        with os.fdopen(fd, "w") as f:
            f.write(content)
        return path

    def test_skips_comment_lines(self):
        path = self._write_temp("# comment\nABC123 desc 123\n")
        rows = parse_phmmer_file(path)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["accession"], "ABC123")
        os.unlink(path)

    def test_skips_blank_lines(self):
        path = self._write_temp("\n\nABC123 desc 50\n\n")
        rows = parse_phmmer_file(path)
        self.assertEqual(len(rows), 1)
        os.unlink(path)

    def test_parses_tlen(self):
        path = self._write_temp("XP_001234 desc 987\n")
        rows = parse_phmmer_file(path)
        self.assertEqual(rows[0]["tlen"], 987)
        os.unlink(path)

    def test_tlen_none_on_non_int(self):
        path = self._write_temp("XP_001234 desc NA\n")
        rows = parse_phmmer_file(path)
        self.assertIsNone(rows[0]["tlen"])
        os.unlink(path)

    def test_skips_short_lines(self):
        path = self._write_temp("onlyone\n")
        rows = parse_phmmer_file(path)
        self.assertEqual(len(rows), 0)
        os.unlink(path)

    def test_line_number_recorded(self):
        path = self._write_temp("# comment\nABC 1 2\nDEF 3 4\n")
        rows = parse_phmmer_file(path)
        self.assertEqual(rows[0]["line_number"], 2)
        self.assertEqual(rows[1]["line_number"], 3)
        os.unlink(path)


class TestPrintProgress(unittest.TestCase):
    def test_returns_incremented_step(self):
        start = time.time()
        # At 50/100, step=10 → should fire at 50%
        result = print_progress("TEST", 50, 100, start, 50, 10)
        self.assertEqual(result, 60)

    def test_no_fire_below_step(self):
        start = time.time()
        result = print_progress("TEST", 5, 100, start, 10, 10)
        self.assertEqual(result, 10)  # step unchanged

    def test_handles_zero_total(self):
        start = time.time()
        result = print_progress("TEST", 0, 0, start, 10, 10)
        self.assertEqual(result, 10)


class TestResolveCacheOut(unittest.TestCase):
    def test_returns_given_path(self):
        self.assertEqual(resolve_cache_out("/tmp/my_cache.json"), "/tmp/my_cache.json")

    def test_generates_path_when_none(self):
        path = resolve_cache_out(None)
        self.assertIn("extraction_cache_", path)
        self.assertTrue(path.endswith(".json"))


if __name__ == "__main__":
    unittest.main()
