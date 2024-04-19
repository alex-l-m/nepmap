'''Look up tridentate carbenes in the Cambridge Structural Database and save the doi of the associated papers'''

import ccdc.search

tridentate_carbene_smarts = 'N~C(~[Ir]12)~N~c:c~1:c~N~C~2~N'

# Search based on the smarts string
substructure = ccdc.search.SMARTSSubstructure(tridentate_carbene_smarts)
search = ccdc.search.SubstructureSearch()
search.add_substructure(substructure)
hits = search.search()

# Print the identifiers of the hits
seen = set()
for hit in hits:
    if hit.identifier not in seen:
        print(hit.identifier)
        seen.add(hit.identifier)
