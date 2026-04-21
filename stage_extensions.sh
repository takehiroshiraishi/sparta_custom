#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CUSTOM_DIR="${ROOT_DIR}/custom/extensions"
SRC_DIR="${ROOT_DIR}/src"
MANIFEST="${ROOT_DIR}/custom/.staged-extension-links"

cleanup_manifest_links() {
  if [[ ! -f "${MANIFEST}" ]]; then
    return
  fi

  while IFS= read -r path; do
    [[ -z "${path}" ]] && continue
    if [[ -L "${path}" ]]; then
      rm -f "${path}"
    fi
  done < "${MANIFEST}"
}

if [[ ! -d "${SRC_DIR}" ]]; then
  echo "src directory not found: ${SRC_DIR}" >&2
  exit 1
fi

cleanup_manifest_links

if [[ ! -d "${CUSTOM_DIR}" ]]; then
  : > "${MANIFEST}"
  exit 0
fi

mapfile -t files < <(
  find "${CUSTOM_DIR}" -type f -path '*/src/*' \
    \( -name '*.cpp' -o -name '*.h' -o -name '*.hpp' \) | sort
)

declare -A basenames=()
staged=()

for source_path in "${files[@]}"; do
  base_name="$(basename "${source_path}")"
  target_path="${SRC_DIR}/${base_name}"

  if [[ -n "${basenames[${base_name}]:-}" ]]; then
    echo "duplicate custom extension filename: ${base_name}" >&2
    echo "  ${basenames[${base_name}]}" >&2
    echo "  ${source_path}" >&2
    exit 1
  fi
  basenames["${base_name}"]="${source_path}"

  if [[ -e "${target_path}" && ! -L "${target_path}" ]]; then
    echo "custom extension conflicts with existing src file: ${target_path}" >&2
    exit 1
  fi

  ln -sfn "${source_path}" "${target_path}"
  staged+=("${target_path}")
done

printf '%s\n' "${staged[@]}" > "${MANIFEST}"

echo "Staged ${#staged[@]} custom extension file(s) into ${SRC_DIR}"
