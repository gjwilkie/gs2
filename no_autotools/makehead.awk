BEGIN {
    print "/* this file automatically created by makehead.awk. */\n"
    ind = 0;
}

{
    if ($1=="type" && $2=="::") {
	print "struct", $3, "{"
	ind += 1;
    } else if (ind > 0) {
	if ($1=="integer" && $2=="::") {
	    for (i=1; i<=ind; i++) {
		printf ("  ");
	    }
	    printf ("int ")
	    for (i=3; i<NF; i++) {
		printf ("%s ", $i);
	    }
	    printf ("%s;\n", $NF);
	} else if ($1=="end" && $2=="type") {
	    print "};\n"
	    ind -= 1;
	}
    }
}
