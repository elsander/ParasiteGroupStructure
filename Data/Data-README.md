## File Names

Each file contains the following information, separated by dashes:
- "Dunne2013-" followed by the food web name
- a flag for concomitant predation (only present in `mat` files)
- the file type

## File Structure

Species in a given food web are ordered in the same way across all file types.

`mat` files files contain the adjacency matrix for the food web,
built from the adjacency lists available from Dunne 2013.

`Category` files contain the ecological categorization of each
species, separated by a newline.

`SpeciesList` files contain the following information,
comma-separated:
- ID (row/column number of the species in the matrix)
- Type (1 for free-living, 2 for parasite)
- Species name
