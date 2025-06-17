import os
import sys
import unittest

# Add src directory to Python path
this_dir = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(this_dir, '..', 'src'))

from pdbstruct.parse import load_soup
from pdbstruct.soup import Soup


class TestParseEquivalence(unittest.TestCase):
    """Test that PDB and CIF parsers produce equivalent soup objects."""

    def setUp(self):
        """Set up test fixtures."""
        self.pdb_file = os.path.join(this_dir, "../examples/1ssx.pdb")
        self.cif_file = os.path.join(this_dir, "../examples/1ssx.cif")

        # Check if test files exist
        if not os.path.exists(self.pdb_file):
            self.skipTest(f"Test file {self.pdb_file} not found")
        if not os.path.exists(self.cif_file):
            self.skipTest(f"Test file {self.cif_file} not found")

    def test_files_exist(self):
        """Test that required example files exist."""
        self.assertTrue(
            os.path.exists(self.pdb_file), f"PDB file {self.pdb_file} should exist"
        )
        self.assertTrue(
            os.path.exists(self.cif_file), f"CIF file {self.cif_file} should exist"
        )

    def test_basic_loading(self):
        """Test that both files can be loaded without errors."""
        pdb_soup = load_soup(self.pdb_file, scrub=True)
        cif_soup = load_soup(self.cif_file, scrub=True)

        self.assertIsInstance(pdb_soup, Soup)
        self.assertIsInstance(cif_soup, Soup)
        self.assertFalse(pdb_soup.is_empty())
        self.assertFalse(cif_soup.is_empty())

    def test_atom_count_equivalence(self):
        """Test that PDB and CIF files have the same number of atoms."""
        pdb_soup = load_soup(self.pdb_file, scrub=True)
        cif_soup = load_soup(self.cif_file, scrub=True)

        pdb_atom_count = pdb_soup.get_atom_count()
        cif_atom_count = cif_soup.get_atom_count()

        self.assertEqual(
            pdb_atom_count,
            cif_atom_count,
            f"Atom counts differ: PDB={pdb_atom_count}, CIF={cif_atom_count}",
        )

    def test_residue_count_equivalence(self):
        """Test that PDB and CIF files have the same number of residues."""
        pdb_soup = load_soup(self.pdb_file, scrub=True)
        cif_soup = load_soup(self.cif_file, scrub=True)

        pdb_res_count = pdb_soup.get_residue_count()
        cif_res_count = cif_soup.get_residue_count()

        self.assertEqual(
            pdb_res_count,
            cif_res_count,
            f"Residue counts differ: PDB={pdb_res_count}, CIF={cif_res_count}",
        )

    def test_atom_coordinates_equivalence(self):
        """Test that atom coordinates are equivalent between PDB and CIF."""
        pdb_soup = load_soup(self.pdb_file, scrub=True)
        cif_soup = load_soup(self.cif_file, scrub=True)

        pdb_atom = pdb_soup.get_atom_proxy()
        cif_atom = cif_soup.get_atom_proxy()

        atom_count = min(pdb_soup.get_atom_count(), cif_soup.get_atom_count())

        for i_atom in range(atom_count):
            pdb_atom.load(i_atom)
            cif_atom.load(i_atom)

            # Compare coordinates with tolerance for floating point precision
            self.assertAlmostEqual(
                pdb_atom.pos.x,
                cif_atom.pos.x,
                places=3,
                msg=f"X coordinate differs for atom {i_atom}",
            )
            self.assertAlmostEqual(
                pdb_atom.pos.y,
                cif_atom.pos.y,
                places=3,
                msg=f"Y coordinate differs for atom {i_atom}",
            )
            self.assertAlmostEqual(
                pdb_atom.pos.z,
                cif_atom.pos.z,
                places=3,
                msg=f"Z coordinate differs for atom {i_atom}",
            )

    def test_atom_properties_equivalence(self):
        """Test that atom properties (B-factor, occupancy, element) are equivalent."""
        pdb_soup = load_soup(self.pdb_file, scrub=True)
        cif_soup = load_soup(self.cif_file, scrub=True)

        pdb_atom = pdb_soup.get_atom_proxy()
        cif_atom = cif_soup.get_atom_proxy()

        atom_count = min(pdb_soup.get_atom_count(), cif_soup.get_atom_count())

        for i_atom in range(atom_count):
            pdb_atom.load(i_atom)
            cif_atom.load(i_atom)

            # Compare B-factors
            self.assertAlmostEqual(
                pdb_atom.bfactor,
                cif_atom.bfactor,
                places=2,
                msg=f"B-factor differs for atom {i_atom}",
            )

            # Compare occupancy
            self.assertAlmostEqual(
                pdb_atom.occupancy,
                cif_atom.occupancy,
                places=2,
                msg=f"Occupancy differs for atom {i_atom}",
            )

            # Compare element
            self.assertEqual(
                pdb_atom.elem,
                cif_atom.elem,
                f"Element differs for atom {i_atom}: PDB={pdb_atom.elem}, CIF={cif_atom.elem}",
            )

            # Compare atom type (may need normalization)
            pdb_atom_type = pdb_atom.atom_type.strip()
            cif_atom_type = cif_atom.atom_type.strip()
            self.assertEqual(
                pdb_atom_type,
                cif_atom_type,
                f"Atom type differs for atom {i_atom}: PDB='{pdb_atom_type}', CIF='{cif_atom_type}'",
            )

    def test_residue_properties_equivalence(self):
        """Test that residue properties are equivalent between PDB and CIF."""
        pdb_soup = load_soup(self.pdb_file, scrub=True)
        cif_soup = load_soup(self.cif_file, scrub=True)

        pdb_res = pdb_soup.get_residue_proxy()
        cif_res = cif_soup.get_residue_proxy()

        res_count = min(pdb_soup.get_residue_count(), cif_soup.get_residue_count())

        for i_res in range(res_count):
            pdb_res.load(i_res)
            cif_res.load(i_res)

            # Compare residue type
            self.assertEqual(
                pdb_res.res_type,
                cif_res.res_type,
                f"Residue type differs for residue {i_res}: PDB={pdb_res.res_type}, CIF={cif_res.res_type}",
            )

            # Compare insertion code
            self.assertEqual(
                pdb_res.ins_code,
                cif_res.ins_code,
                f"Insertion code differs for residue {i_res}: PDB={pdb_res.ins_code}, CIF={cif_res.ins_code}",
            )

    def test_structure_metadata_equivalence(self):
        """Test that structure metadata is equivalent."""
        pdb_soup = load_soup(self.pdb_file, scrub=True)
        cif_soup = load_soup(self.cif_file, scrub=True)

        # Both should have structure IDs
        self.assertTrue(
            len(pdb_soup.structure_ids) > 0, "PDB soup should have structure IDs"
        )
        self.assertTrue(
            len(cif_soup.structure_ids) > 0, "CIF soup should have structure IDs"
        )

        # Structure ID should be based on filename (1ssx)
        pdb_base_id = pdb_soup.structure_ids[0].split("[")[
            0
        ]  # Remove model numbers if present
        cif_base_id = cif_soup.structure_ids[0].split("[")[0]

        self.assertEqual(
            pdb_base_id,
            cif_base_id,
            f"Base structure IDs differ: PDB={pdb_base_id}, CIF={cif_base_id}",
        )

    def test_atom_residue_mapping_equivalence(self):
        """Test that atom-to-residue mapping is equivalent."""
        pdb_soup = load_soup(self.pdb_file, scrub=True)
        cif_soup = load_soup(self.cif_file, scrub=True)

        pdb_atom = pdb_soup.get_atom_proxy()
        cif_atom = cif_soup.get_atom_proxy()

        atom_count = min(pdb_soup.get_atom_count(), cif_soup.get_atom_count())

        for i_atom in range(atom_count):
            pdb_atom.load(i_atom)
            cif_atom.load(i_atom)

            # The residue index might differ, but the residue properties should match
            pdb_res = pdb_soup.get_residue_proxy(pdb_atom.i_res)
            cif_res = cif_soup.get_residue_proxy(cif_atom.i_res)

            self.assertEqual(
                pdb_res.res_type,
                cif_res.res_type,
                f"Mapped residue type differs for atom {i_atom}: {i_atom}: {pdb_atom} <-> {cif_atom}",
            )

    def test_scrub_option_equivalence(self):
        """Test that scrub option produces equivalent results for both formats."""
        pdb_soup_scrub = load_soup(self.pdb_file, scrub=True)
        cif_soup_scrub = load_soup(self.cif_file, scrub=True)
        pdb_soup_no_scrub = load_soup(self.pdb_file, scrub=False)
        cif_soup_no_scrub = load_soup(self.cif_file, scrub=False)

        # Scrubbed versions should have same or fewer atoms
        self.assertLessEqual(
            pdb_soup_scrub.get_atom_count(), pdb_soup_no_scrub.get_atom_count()
        )
        self.assertLessEqual(
            cif_soup_scrub.get_atom_count(), cif_soup_no_scrub.get_atom_count()
        )

        # Scrubbed versions should be equivalent to each other
        self.assertEqual(
            pdb_soup_scrub.get_atom_count(), cif_soup_scrub.get_atom_count()
        )
        self.assertEqual(
            pdb_soup_scrub.get_residue_count(), cif_soup_scrub.get_residue_count()
        )


if __name__ == "__main__":
    unittest.main()
