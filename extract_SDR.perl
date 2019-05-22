#!/usr/bin/perl

open FH, "<SN_M_10TO150_count_tab_025cov_pos_diversity.depth";
open OUT1, ">SN_M_10TO150_count_tab_025cov_pos_diversity_SNsdr.depth";
open OUT2, ">SN_M_10TO150_count_tab_025cov_pos_diversity_SEsdr.depth";
$n = 0;

while (<FH>)
{
    $seq = $_;
    chomp $seq;
    $n = $n + 1;
    if ($n == 1)
    {
	print OUT1 "$seq\n";
	print OUT2 "$seq\n";
    }
    elsif (/^\S+\s+(\S+)\s+(\S+)/)
    {
	$chr = $1;
	$pos = $2;
	if ($chr eq "Chr07")
	{
	    if (($pos >= 4886044) && ($pos <= 6857240))
	    {
		print OUT1 "$seq\n";
	    }
	}
	elsif ($chr eq "Chr15W")
	{
	    if (($pos >= 6131694) && ($pos <= 8359477))
	    {
		print OUT2 "$seq\n";
	    }
	}
    }
}
close FH;
close OUT1;
close OUT2;
