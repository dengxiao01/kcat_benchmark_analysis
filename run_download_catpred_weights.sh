#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

CATPRED_CAPSULE_DIR="${CATPRED_CAPSULE_DIR:-$ROOT/external_methods/CatPred_capsule}"
ARCHIVE="$CATPRED_CAPSULE_DIR/capsule_data_update.tar.gz"
PART="$CATPRED_CAPSULE_DIR/capsule_data_update.tar.gz.part"
EXPECTED_BYTES=10207582747

URL_PRIMARY="https://catpred.s3.us-east-1.amazonaws.com/capsule_data_update.tar.gz"
URL_FALLBACK="https://catpred.s3.amazonaws.com/capsule_data_update.tar.gz"

mkdir -p "$CATPRED_CAPSULE_DIR"

download_one() {
  local url="$1"
  curl \
    --location \
    --fail \
    --continue-at - \
    --retry 20 \
    --retry-all-errors \
    --retry-delay 10 \
    --connect-timeout 30 \
    --speed-time 120 \
    --speed-limit 10240 \
    --output "$PART" \
    "$url"
}

size_of() {
  if [[ -f "$1" ]]; then
    stat -c '%s' "$1"
  else
    echo 0
  fi
}

if [[ -f "$ARCHIVE" ]] && [[ "$(size_of "$ARCHIVE")" -eq "$EXPECTED_BYTES" ]]; then
  echo "Archive already complete: $ARCHIVE"
else
  echo "Downloading CatPred capsule to: $PART"
  echo "Expected bytes: $EXPECTED_BYTES"
  download_one "$URL_PRIMARY" || download_one "$URL_FALLBACK"

  actual_bytes="$(size_of "$PART")"
  if [[ "$actual_bytes" -ne "$EXPECTED_BYTES" ]]; then
    echo "Download is incomplete: $actual_bytes bytes, expected $EXPECTED_BYTES bytes." >&2
    echo "Re-run this script to resume." >&2
    exit 2
  fi

  mv -f "$PART" "$ARCHIVE"
fi

echo "Checking archive integrity..."
tar -tzf "$ARCHIVE" >/dev/null

if [[ "${EXTRACT:-1}" != "0" ]]; then
  echo "Extracting into: $CATPRED_CAPSULE_DIR"
  tar -xzf "$ARCHIVE" -C "$CATPRED_CAPSULE_DIR"
fi

KCAT_DIR="$CATPRED_CAPSULE_DIR/data/pretrained/production/kcat"
if [[ -d "$KCAT_DIR" ]]; then
  echo "CatPred kcat checkpoint ready: $KCAT_DIR"
else
  echo "Archive extracted, but kcat checkpoint directory was not found at: $KCAT_DIR" >&2
  echo "Inspect the extracted data directory before running CatPred." >&2
  exit 2
fi
