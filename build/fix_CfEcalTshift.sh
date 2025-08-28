#!/bin/bash
set -euo pipefail
FILE="CfEcalTshift.scr"
[ -f "$FILE" ] || { echo "Fichier non trouvé: $FILE"; exit 1; }

cp -a -- "$FILE" "${FILE}.bak"

awk '
BEGIN {
  inblk=0; bufn=0; indent=""
}

# début du tableau runlists=( ... )
/^[[:space:]]*runlists[[:space:]]*=\s*\(/ {
  print
  inblk=1
  indent=""
  next
}

# fin du tableau
inblk && /^[[:space:]]*\)/ {
  for (i=1; i<=bufn; i++) {
    raw = buf[i]

    # décommente le début de ligne et trim gauche
    line = raw
    sub(/^[[:space:]]*#+[[:space:]]*/, "", line)
    sub(/^[[:space:]]*/, "", line)

    # extrait le nom entre guillemets pour la clé de déduplication
    key=""
    if (match(line, /"([^"]+)"/, m)) {
      key = m[1]
    }

    if (key ~ /^run[0-9]+EcalTshift\.dat$/) {
      if (!seen[key]++) {
        print indent line
      }
    } else {
      if (line ~ /[^[:space:]]/) {
        print indent line
      } else {
        print
      }
    }
  }
  print
  inblk=0; bufn=0
  next
}

# à l intérieur du bloc, on stocke les lignes (et on capture l indent au 1er passage)
inblk {
  if (bufn==0) {
    if (match($0, /^[[:space:]]+/, im)) indent=im[0]; else indent=""
  }
  buf[++bufn] = $0
  next
}

# en dehors du bloc, impression normale
{ print }
' "$FILE" > "${FILE}.tmp" && mv -- "${FILE}.tmp" "$FILE"

echo "Fini. Sauvegarde: ${FILE}.bak"