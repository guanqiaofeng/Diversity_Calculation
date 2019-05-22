#!/usr/bin/perl

$n = 0;
$pi00 = ();
$theta00 = ();

open FH, "<SN_M_10TO150_count_tab_025cov_pos_diversity.depth";
open OUT, ">SN_M_10TO150_count_tab_025cov_pos_diversity_persite.depth";
print OUT "theta_per_site\tpi_per_site\n";
while (<FH>)
{
    $n = $n + 1;
    if (/^\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)/)
    {
	$theta = $1;
	$pi = $2;
	if ($n == 1000)
	{
	    $theta00 = $theta00 / $n;
	    $pi00 = $pi00 / $n;
	    print OUT "$theta00\t$pi00\n";
	    $n = 0;
	    $theta00 = $theta;
	    $pi00 = $pi;
	}
	else
	{
	    $theta00 = $theta00 + $theta;
	    $pi00 = $pi00 + $pi;
	}
    }
}
$theta00 = $theta00 / $n;
$pi00 = $pi00 / $n;
print OUT "$theta00\t$pi00\n";
close FH;
close OUT;
