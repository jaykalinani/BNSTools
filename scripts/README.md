# StellarCollapse Table workflow
1. Download table from https://stellarcollapse.org/microphysics.html
2. If desired, extend the table using `extend_table_V3.py`. Example usage:
   ```
   python extend_table_V3.py -t <your table>.h5 --min-dens -3 --stitch-dens 3.4 --width 0.75 --tol 3e-2 --min-logtemp 999 --stitch-temp 999
   ```
   See script for descriptions of arguments. The new table will be saved as `<your_table>_ext.h5`.
3. Adjust table using `python adjust_table_mineps0.py <your table>.h5`
   - This fixes unphysical, NaN, or infinite sound speeds.
   - This also transfers energy from rest baryon mass to internal energy such that the minimum specific internal energy (epsilon) is zero. Record the new baryon mass, as it is critical to keep this value consistent when generating and loading initial data.
   - The adjusted table is saved as `<your table>_adj_mineps0.h5`.
4. If desired, soften the high density end of the table for triggered collapse using `soften_table.py`. Example usage:
   ```
   python soften_table.py -t <your table>.h5 -s <soft_table>.h5 --stitch-dens 14.65 --width 0.2 --tol 3e-2
   ```
   See script for descriptions of arguments. The new table will be saved as `<your table>_collapse.h5`.
5. Use the EOS tools from the THC repository (https://bitbucket.org/FreeTHC/thcextra/src/master/tools/eos_tables/) to generate 1D slices of the table for initial data. Change the baryon mass defined in `modules/constants.py` to your adjusted value.
6. Use `python generate_Ye_of_rho.py <your table>.h5 <temperature in MeV>` to generate an auxillary table of `Ye` values in beta equilibrium used by `ID_TabEOS_HydroQuantities`. This script currently only supports constant temperature EOS slices. 
