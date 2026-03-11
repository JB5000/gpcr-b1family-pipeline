#!/usr/bin/env python3
"""
Unit tests for build_excels_from_cache.py
Covers the pure helper functions:
  - sanitize_name
  - normalize_for_match
  - classify_subfamily
  - get_species_abbrev
  - is_np
  - normalize_sheet_title
"""

import unittest

from build_excels_from_cache import (
    classify_subfamily,
    get_species_abbrev,
    is_np,
    normalize_for_match,
    normalize_sheet_title,
    sanitize_name,
)


class TestSanitizeName(unittest.TestCase):
    def test_alphanumeric_preserved(self):
        self.assertEqual(sanitize_name("GCGR"), "GCGR")

    def test_spaces_replaced(self):
        result = sanitize_name("my file name")
        self.assertNotIn(" ", result)

    def test_hyphens_replaced(self):
        result = sanitize_name("GLP-1R")
        self.assertNotIn("-", result)

    def test_empty_fallback(self):
        self.assertEqual(sanitize_name(""), "output")

    def test_only_symbols_fallback(self):
        self.assertEqual(sanitize_name("@#$"), "output")


class TestNormalizeForMatch(unittest.TestCase):
    def test_lowercases(self):
        self.assertEqual(normalize_for_match("GCGR"), "gcgr")

    def test_strips_non_alphanumeric(self):
        result = normalize_for_match("glucagon-like peptide-1")
        self.assertNotIn("-", result)

    def test_collapses_spaces(self):
        result = normalize_for_match("  a  b  ")
        self.assertEqual(result, "a b")


class TestClassifySubfamily(unittest.TestCase):
    def test_gcgr_keyword(self):
        desc = "glucagon receptor [Homo sapiens]"
        self.assertEqual(classify_subfamily(desc), "GCGR")

    def test_glp1_keyword(self):
        desc = "glucagon-like peptide-1 receptor [Mus musculus]"
        self.assertEqual(classify_subfamily(desc), "GLP1R")

    def test_crhr1_keyword(self):
        desc = "corticotropin releasing hormone receptor 1 [Danio rerio]"
        self.assertEqual(classify_subfamily(desc), "CRHR1")

    def test_pdfr_accession_override(self):
        # This accession is explicitly mapped to PDFR regardless of description
        self.assertEqual(classify_subfamily("unrelated desc", "XP_052822567"), "PDFR")

    def test_pdfr_accession_with_version(self):
        self.assertEqual(classify_subfamily("unrelated desc", "XP_052822567.1"), "PDFR")

    def test_unclassified_unknown(self):
        desc = "hypothetical protein XP_99999 [Unknown species]"
        self.assertEqual(classify_subfamily(desc), "UNCLASSIFIED")

    def test_sctr_keyword(self):
        desc = "secretin receptor [Gallus gallus]"
        self.assertEqual(classify_subfamily(desc), "SCTR")

    def test_pth1r_keyword(self):
        desc = "parathyroid hormone receptor 1 [Xenopus laevis]"
        self.assertEqual(classify_subfamily(desc), "PTH1R")


class TestGetSpeciesAbbrev(unittest.TestCase):
    def test_two_word_species(self):
        desc = "secretin receptor [Homo sapiens]"
        self.assertEqual(get_species_abbrev(desc), "Hsa")

    def test_single_word_species(self):
        desc = "receptor [Drosophila]"
        self.assertEqual(get_species_abbrev(desc), "NA")

    def test_no_brackets(self):
        desc = "receptor without species"
        self.assertEqual(get_species_abbrev(desc), "NA")

    def test_last_bracket_used(self):
        desc = "receptor [old] final desc [Danio rerio]"
        abbrev = get_species_abbrev(desc)
        self.assertTrue(abbrev.startswith("D"))

    def test_multiword_uses_first_and_last(self):
        desc = "receptor [Takifugu rubripes]"
        self.assertEqual(get_species_abbrev(desc), "Tru")


class TestIsNp(unittest.TestCase):
    def test_np_accession(self):
        self.assertTrue(is_np("NP_001234"))

    def test_xp_accession(self):
        self.assertFalse(is_np("XP_001234"))

    def test_empty_string(self):
        self.assertFalse(is_np(""))


class TestNormalizeSheetTitle(unittest.TestCase):
    def test_truncates_long_title(self):
        long_title = "A" * 40
        used: set[str] = set()
        result = normalize_sheet_title(long_title, used)
        self.assertLessEqual(len(result), 31)

    def test_deduplicates(self):
        used: set[str] = set()
        t1 = normalize_sheet_title("SameSheet", used)
        t2 = normalize_sheet_title("SameSheet", used)
        self.assertNotEqual(t1, t2)

    def test_first_use_unchanged_if_short(self):
        used: set[str] = set()
        result = normalize_sheet_title("ShortTitle", used)
        self.assertEqual(result, "ShortTitle")


if __name__ == "__main__":
    unittest.main()
