# ...existing code...
#!/bin/bash
# Crée Todo.dat contenant la liste "run414.dat" .. "run536.dat", une par ligne

OUT="/data/Malia/Analyse_FASTER/build/Todo.dat"

# vide/cree le fichier de sortie
: > "$OUT"

for n in $(seq 414 536); do
  printf '"run%d.dat\n"' "$n" >> "$OUT"
done

echo "Todo créé: $OUT"
# ...existing code...