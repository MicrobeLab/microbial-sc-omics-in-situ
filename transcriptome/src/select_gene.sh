# $1 polyA position file
# $2 contig id
# $3 gene start postion
# $4 gene end position
# $5 gene annotation

awk -F "," 'BEGIN{before=0;after=0;print "Gene annotation\tid\tgene_start\tgene_end\tpolyA_1_2_3...\tbefore_pos\tafter_pos\n"}{

if( $1=="'$2'" ){  printf("%s","'$5'\t'$2'\t'$3'\t'$4'"); 
for(t=3;t<=NF;t++){if($t< '$4') { if($t> '$3' ) {printf("%s","\t"$t )  } }
		if( $t < '$3') { before=$t}
			if( $t > '$4'){ if(after==0) {after=$t;}}	
		}
		printf("%s","\t"before"\t"after"\n")
		}
	}' $1
