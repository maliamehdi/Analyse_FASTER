#!/bin/bash
# Supprime les 3 dernières lignes de run414.dat .. run672.dat
# Usage: chmod +x trim_runs.sh && ./trim_runs.sh
DIR="."            # chemin vers le dossier contenant les .dat (modifier si besoin)
MAKE_BACKUP=0      # 1 pour créer fichier.bak avant modification, 0 sinon

for n in $(seq 536 672); do
  f="${DIR}/run${n}.dat"
  [ -f "$f" ] || continue
  if [ "$MAKE_BACKUP" -eq 1 ]; then
    cp -a -- "$f" "${f}.bak"
  fi
  # head -n -3 fonctionne sur GNU head (Linux). Remplace le fichier par son contenu sans les 3 dernières lignes.
  head -n -2 -- "$f" > "${f}.tmp" && mv -- "${f}.tmp" "$f"
  echo "trimmed: $f"
done