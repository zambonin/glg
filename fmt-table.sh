#!/usr/bin/env sh

function pprint() {
  mlr \
    --opprint \
    clean-whitespace \
    then cut -f p,d,n,alg,"$1" \
    then reshape -s alg,"$1" \
    then put '
      var min_val = "";
      for (k, v in $*) {
        if (k != "p" && k != "d" && k != "n" && is_numeric(v)) {
          if (min_val == "") {
            min_val = v;
          } else {
            min_val = min(min_val, v);
          }
        }
      }

      if (min_val != "") {
        var new_record = {"p": $p, "d": $d, "n": $n};

        for (k, v in $*) {
          if (k != "p" && k != "d" && k != "n" && is_numeric(v)) {
            new_record[k] = v;
            new_record[k . "_ratio"]
              = (v == min_val ? "*" : fmtnum(v / min_val, "%.2f"));
          }
        }
        $* = new_record;
      }
    ' "$2"
}

function pprint_cycles() {
  pprint avg_cycles "$1"
}

function pprint_limbs() {
  pprint limb "$1"
}

INPUT="/dev/stdin"
OPTION=

while getopts clf:-: OPT ; do
  if [ "$OPT" = "-" ] ; then
    OPT="${OPTARG%%=*}"
    OPTARG="${OPTARG#$OPT}"
    OPTARG="${OPTARG#=}"
  fi
  case "$OPT" in
    c | cycles ) OPTION="cycles" ;;
    f | file ) INPUT="${OPTARG:-$INPUT}" ;;
    l | limbs ) OPTION="limbs" ;;
    ??* ) exit 1 ;;
    \? ) exit 2 ;;
  esac
done
shift $((OPTIND - 1))

if [ -z "$OPTION" ] || [ ! -r "$INPUT" ] ; then
  cat <<-EOF
Usage: sh ${0##*/} [options]
  -c,  --cycles
  -f,  --file=<path>
  -l,  --limbs
EOF
  exit 1
fi

case "$OPTION" in
  cycles ) pprint_cycles "$INPUT" ;;
  limbs ) pprint_limbs "$INPUT" ;;
esac
