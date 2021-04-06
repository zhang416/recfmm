#!/usr/bin/env bash

EXEC=$1; shift

$EXEC "$@" | grep "\*\{7\} pass \*\{7\}"
