cat  $1 | awk '{for (i=1; i <= 6; i++) printf("%s ",$i); printf(" @ 1.2 16 C [red]\n");}'
