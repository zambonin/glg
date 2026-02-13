#!/usr/bin/env sh

function pprint_cycles() {
  mlr \
    --opprint \
    clean-whitespace \
    then put '$q = $p**$d' \
    then cut -o -f n,q,alg,avg_cycles \
    then sort -n n \
    then reshape -s alg,avg_cycles /dev/stdin
}

function pprint_limbs() {
  mlr \
    --opprint \
    clean-whitespace \
    then put '$q = $p**$d' \
    then cut -o -f n,q,chunk,alg,limb \
    then sort -n n \
    then filter '$chunk == 1 || $chunk == 4 || $chunk == 16 || $chunk == 64' \
    then filter '$n == 5 || $n == 11 || $n == 30' \
    then reshape -s alg,limb /dev/stdin
}

OPTION=

while getopts clf:-: OPT ; do
  if [ "$OPT" = "-" ] ; then
    OPT="${OPTARG%%=*}"
    OPTARG="${OPTARG#$OPT}"
    OPTARG="${OPTARG#=}"
  fi
  case "$OPT" in
    c | cycles ) OPTION="cycles" ;;
    l | limbs ) OPTION="limbs" ;;
    ??* ) exit 1 ;;
    \? ) exit 2 ;;
  esac
done
shift $((OPTIND - 1))

if [ -z "$OPTION" ] ; then
  cat <<-EOF
Usage: sh ${0##*/} [options]
  -c,  --cycles
  -l,  --limbs
EOF
  exit 1
fi

case "$OPTION" in
  cycles )
    make --silent clean prof-c | pprint_cycles ;;
  limbs )
    for i in 1 4 16 64 ; do
      make CHUNK_SIZE="$i" --silent clean prof-c
    done | pprint_limbs ;;
esac
