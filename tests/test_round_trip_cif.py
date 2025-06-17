import os
import tempfile
import unittest
import sys

# Add src directory to Python path
this_dir = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(this_dir, '../src'))

from pdbstruct.parse import load_soup, write_soup
from pdbstruct.soup import Soup


class TestCifRoundtrip(unittest.TestCase):
    """Test CIF file reading, writing, and re-reading consistency."""

    def setUp(self):
        """Set up test fixtures."""
        self.cif_file = os.path.join(this_dir, "../examples/1ssx.cif")
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up temporary files."""
        import shutil

        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_cif_file_exists(self):
        """Test that the example CIF file exists."""
        self.assertTrue(
            os.path.exists(self.cif_file), f"Example CIF file {self.cif_file} not found"
        )

    def compare_soups(self, soup1: Soup, soup2: Soup, test_name: str):
        """Compare two soup objects for equality."""
        # Basic counts
        self.assertEqual(
            soup1.get_atom_count(),
            soup2.get_atom_count(),
            f"{test_name}: Atom counts differ",
        )
        self.assertEqual(
            soup1.get_residue_count(),
            soup2.get_residue_count(),
            f"{test_name}: Residue counts differ",
        )

        # Structure metadata
        self.assertEqual(
            len(soup1.structure_ids),
            len(soup2.structure_ids),
            f"{test_name}: Structure ID counts differ",
        )
        self.assertEqual(
            len(soup1.chains), len(soup2.chains), f"{test_name}: Chain counts differ"
        )

        # Value tables
        self.assertEqual(
            len(soup1.elems), len(soup2.elems), f"{test_name}: Element counts differ"
        )
        self.assertEqual(
            len(soup1.atom_types),
            len(soup2.atom_types),
            f"{test_name}: Atom type counts differ",
        )
        self.assertEqual(
            len(soup1.res_types),
            len(soup2.res_types),
            f"{test_name}: Residue type counts differ",
        )

        # Compare atoms
        atom1 = soup1.get_atom_proxy()
        atom2 = soup2.get_atom_proxy()

        for i_atom in range(soup1.get_atom_count()):
            atom1.load(i_atom)
            atom2.load(i_atom)

            # Position comparison with tolerance
            self.assertAlmostEqual(
                atom1.pos.x,
                atom2.pos.x,
                places=3,
                msg=f"{test_name}: Atom {i_atom} X coordinate differs",
            )
            self.assertAlmostEqual(
                atom1.pos.y,
                atom2.pos.y,
                places=3,
                msg=f"{test_name}: Atom {i_atom} Y coordinate differs",
            )
            self.assertAlmostEqual(
                atom1.pos.z,
                atom2.pos.z,
                places=3,
                msg=f"{test_name}: Atom {i_atom} Z coordinate differs",
            )

            # Other atom properties
            self.assertAlmostEqual(
                atom1.bfactor,
                atom2.bfactor,
                places=2,
                msg=f"{test_name}: Atom {i_atom} B-factor differs",
            )
            self.assertAlmostEqual(
                atom1.occupancy,
                atom2.occupancy,
                places=2,
                msg=f"{test_name}: Atom {i_atom} occupancy differs",
            )
            self.assertEqual(
                atom1.elem, atom2.elem, f"{test_name}: Atom {i_atom} element differs"
            )
            self.assertEqual(
                atom1.atom_type,
                atom2.atom_type,
                f"{test_name}: Atom {i_atom} type differs",
            )
            self.assertEqual(
                atom1.alt, atom2.alt, f"{test_name}: Atom {i_atom} alt location differs"
            )
            self.assertEqual(
                atom1.i_res,
                atom2.i_res,
                f"{test_name}: Atom {i_atom} residue index differs",
            )

        # Compare residues
        res1 = soup1.get_residue_proxy()
        res2 = soup2.get_residue_proxy()

        for i_res in range(soup1.get_residue_count()):
            res1.load(i_res)
            res2.load(i_res)

            self.assertEqual(
                res1.res_type,
                res2.res_type,
                f"{test_name}: Residue {i_res} type differs",
            )
            self.assertEqual(
                res1.chain, res2.chain, f"{test_name}: Residue {i_res} chain differs"
            )
            self.assertEqual(
                res1.res_num,
                res2.res_num,
                f"{test_name}: Residue {i_res} number differs",
            )
            self.assertEqual(
                res1.ins_code,
                res2.ins_code,
                f"{test_name}: Residue {i_res} insertion code differs",
            )
            self.assertEqual(
                res1.is_polymer,
                res2.is_polymer,
                f"{test_name}: Residue {i_res} polymer flag differs",
            )

    def test_cif_roundtrip_cif_format(self):
        """Test CIF -> CIF -> CIF roundtrip."""
        # Load original CIF
        original_soup = load_soup(self.cif_file)

        # Write to CIF format
        temp_cif1 = os.path.join(self.temp_dir, "temp1.cif")
        write_soup(original_soup, temp_cif1)

        # Read back the written CIF
        reloaded_soup = load_soup(temp_cif1)

        # Compare soups
        self.compare_soups(original_soup, reloaded_soup, "CIF->CIF roundtrip")

        # Write again to test consistency
        temp_cif2 = os.path.join(self.temp_dir, "temp2.cif")
        write_soup(reloaded_soup, temp_cif2)

        # Read the second written file
        reloaded_soup2 = load_soup(temp_cif2)

        # Compare second roundtrip
        self.compare_soups(reloaded_soup, reloaded_soup2, "CIF->CIF->CIF roundtrip")

    def test_cif_to_pdb_roundtrip(self):
        """Test CIF -> PDB -> PDB roundtrip."""
        # Load original CIF
        original_soup = load_soup(self.cif_file)

        # Write to PDB format
        temp_pdb1 = os.path.join(self.temp_dir, "temp1.pdb")
        write_soup(original_soup, temp_pdb1)

        # Read back the written PDB
        reloaded_soup = load_soup(temp_pdb1)

        # Compare soups (PDB format may have some precision loss)
        self.compare_soups(original_soup, reloaded_soup, "CIF->PDB roundtrip")

        # Write PDB again
        temp_pdb2 = os.path.join(self.temp_dir, "temp2.pdb")
        write_soup(reloaded_soup, temp_pdb2)

        # Read the second written file
        reloaded_soup2 = load_soup(temp_pdb2)

        # Compare second roundtrip
        self.compare_soups(reloaded_soup, reloaded_soup2, "PDB->PDB roundtrip")

    def test_atom_data_integrity(self):
        """Test specific atom data integrity during roundtrip."""
        original_soup = load_soup(self.cif_file)

        # Write and reload
        temp_file = os.path.join(self.temp_dir, "integrity_test.cif")
        write_soup(original_soup, temp_file)
        reloaded_soup = load_soup(temp_file)

        # Test specific atom properties
        if original_soup.get_atom_count() > 0:
            atom_orig = original_soup.get_atom_proxy(0)
            atom_reload = reloaded_soup.get_atom_proxy(0)

            # Test first atom in detail
            self.assertAlmostEqual(atom_orig.pos.x, atom_reload.pos.x, places=3)
            self.assertAlmostEqual(atom_orig.pos.y, atom_reload.pos.y, places=3)
            self.assertAlmostEqual(atom_orig.pos.z, atom_reload.pos.z, places=3)
            self.assertEqual(atom_orig.elem, atom_reload.elem)
            self.assertEqual(atom_orig.atom_type, atom_reload.atom_type)

            # Test last atom
            if original_soup.get_atom_count() > 1:
                last_idx = original_soup.get_atom_count() - 1
                atom_orig.load(last_idx)
                atom_reload.load(last_idx)

                self.assertAlmostEqual(atom_orig.pos.x, atom_reload.pos.x, places=3)
                self.assertAlmostEqual(atom_orig.pos.y, atom_reload.pos.y, places=3)
                self.assertAlmostEqual(atom_orig.pos.z, atom_reload.pos.z, places=3)

    def test_residue_data_integrity(self):
        """Test specific residue data integrity during roundtrip."""
        original_soup = load_soup(self.cif_file)

        # Write and reload
        temp_file = os.path.join(self.temp_dir, "residue_test.cif")
        write_soup(original_soup, temp_file)
        reloaded_soup = load_soup(temp_file)

        # Test specific residue properties
        if original_soup.get_residue_count() > 0:
            res_orig = original_soup.get_residue_proxy(0)
            res_reload = reloaded_soup.get_residue_proxy(0)

            # Test first residue
            self.assertEqual(res_orig.res_type, res_reload.res_type)
            self.assertEqual(res_orig.chain, res_reload.chain)
            self.assertEqual(res_orig.res_num, res_reload.res_num)
            self.assertEqual(res_orig.ins_code, res_reload.ins_code)

            # Test atom counts in residue
            atoms_orig = res_orig.get_atom_indices()
            atoms_reload = res_reload.get_atom_indices()
            self.assertEqual(len(atoms_orig), len(atoms_reload))

    def test_chain_consistency(self):
        """Test chain information consistency during roundtrip."""
        original_soup = load_soup(self.cif_file)

        # Write and reload
        temp_file = os.path.join(self.temp_dir, "chain_test.cif")
        write_soup(original_soup, temp_file)
        reloaded_soup = load_soup(temp_file)

        # Compare chain lists
        self.assertEqual(original_soup.chains, reloaded_soup.chains)

        # Test chain assignments for each residue
        res_orig = original_soup.get_residue_proxy()
        res_reload = reloaded_soup.get_residue_proxy()

        for i_res in range(original_soup.get_residue_count()):
            res_orig.load(i_res)
            res_reload.load(i_res)
            self.assertEqual(
                res_orig.chain, res_reload.chain, f"Chain mismatch at residue {i_res}"
            )

    def test_empty_soup_handling(self):
        """Test handling of empty soup objects."""
        empty_soup = Soup()

        # Write empty soup
        temp_file = os.path.join(self.temp_dir, "empty_test.cif")
        write_soup(empty_soup, temp_file)

        # Read it back
        reloaded_soup = load_soup(temp_file)

        # Should still be empty
        self.assertEqual(reloaded_soup.get_atom_count(), 0)
        self.assertEqual(reloaded_soup.get_residue_count(), 0)

    def test_file_format_detection(self):
        """Test that file format is correctly detected by extension."""
        original_soup = load_soup(self.cif_file)

        # Test .cif extension
        temp_cif = os.path.join(self.temp_dir, "test.cif")
        write_soup(original_soup, temp_cif)

        # Test .pdb extension
        temp_pdb = os.path.join(self.temp_dir, "test.pdb")
        write_soup(original_soup, temp_pdb)

        # Both should be readable
        soup_from_cif = load_soup(temp_cif)
        soup_from_pdb = load_soup(temp_pdb)

        self.assertEqual(soup_from_cif.get_atom_count(), original_soup.get_atom_count())
        self.assertEqual(soup_from_pdb.get_atom_count(), original_soup.get_atom_count())


if __name__ == "__main__":
    # Run tests
    unittest.main()
