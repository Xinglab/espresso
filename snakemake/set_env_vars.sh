#!/bin/bash

function set_conda_env_prefix() {
  local ORIG_DIR="$(pwd)" || return 1

  local REL_SCRIPT_DIR="$(dirname ${BASH_SOURCE[0]})" || return 1
  cd "${REL_SCRIPT_DIR}" || return 1
  local SCRIPT_DIR="$(pwd)" || return 1

  cd "${ORIG_DIR}" || return 1

  CONDA_ENV_PREFIX="${SCRIPT_DIR}/conda_env"
}

function main() {
  # need to use the setup that conda init writes to .bashrc
  source "${HOME}/.bashrc" || return 1
  set_conda_env_prefix || return 1
}

main "$@"
