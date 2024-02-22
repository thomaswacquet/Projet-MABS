

## Utilisation de la librarie
```python
aligner = Aligner("ACGTGCTGG", "ACGTGC", {"identity": 1, "substitution": -1, "indel": -2})
alignment = aligner.FindAlignment()
print(alignment.alignedU)
print(alignment.alignedV)
```

