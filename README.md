# Projet MABS
Théo Demarque
Thomas Wacquet

## Exécution du programme
```
python3 main.py -f data\sequences.fa -c data\config.txt
```

## Utilisation de la librarie
```python
aligner = Aligner("ACGTGCTGG", "ACGTGC", {"identity": 1, "substitution": -1, "indel": -2})
alignment = aligner.FindAlignment()
print(alignment.alignedU)
print(alignment.alignedV)
```