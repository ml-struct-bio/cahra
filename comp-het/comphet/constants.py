"""Ground-truth species in the challenge dataset.

Four target complexes plus two decoy ("junk") species that were deliberately
mixed in and were never meant to be reconstructed as targets. Only the four
targets are scored; the decoys matter because particles of theirs landing in a
target class drag that class's purity down.
"""

TARGET_SPECIES = ["PKM2", "GDH", "ADH", "Aldolase"]
JUNK_SPECIES = ["Gamma-Globulin", "Beta-Amylase"]
SPECIES_ORDER = TARGET_SPECIES + JUNK_SPECIES

# The raw source.star labels the gamma globulin sample "IgG", which is not what
# it is -- the sample is gamma globulin. Renamed on load so the correct name is
# what appears in every figure and CSV.
SPECIES_RENAME = {"IgG": "Gamma-Globulin"}
