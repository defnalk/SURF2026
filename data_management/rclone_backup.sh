#!/usr/bin/env bash
# rclone_backup.sh -- SURF 2026
# Nightly one-way sync of local project directory to Purdue Anvil.
# Scheduled by com.defne.surf.backup.plist (launchd).

set -euo pipefail

LOCAL_DIR="$HOME/Desktop/defnalk/SURF2026"
REMOTE="anvil:/anvil/scratch/x-USER/SURF2026_backup"   # edit remote name
LOG="$HOME/Library/Logs/surf2026_rclone.log"

mkdir -p "$(dirname "$LOG")"

{
    echo "---- $(date) ----"
    rclone sync "$LOCAL_DIR" "$REMOTE" \
        --exclude '.git/**' \
        --exclude '*.dcd' --exclude '*.dump' --exclude '*.lammpstrj' \
        --exclude '__pycache__/**' \
        --transfers 4 --checkers 8 \
        --log-level INFO
    echo "---- done $(date) ----"
} >> "$LOG" 2>&1
