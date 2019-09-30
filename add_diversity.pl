#!/usr/bin/perl

$a18 = 4.14678;  # 1 + 1/2 + ... + 1/35
$a19 = 4.20159;  # 1 + 1/2 + ... + 1/37
$a20 = 4.25354;  # 1 + 1/2 + ... + 1/39
$a21 = 4.30293;  # 1 + 1/2 + ... + 1/41
$a22 = 4.35000;  # 1 + 1/2 + ... + 1/43
$a23 = 4.39495;  # 1 + 1/2 + ... + 1/45
$a24 = 4.43796;  # 1 + 1/2 + ... + 1/47 

%a = ();

open FH, "<SN_M_10TO150_count_tab_025cov.depth";
while (<FH>)
{
    if (/(\S+)\s+(\S+)\s+(\S+)/)
    {
	$cov = $1;
	$chr = $2;
	$pos = $3;
	$chrpos = $chr . "_" . $pos;
	if ($cov == 18)
	{
	    $a{$chrpos} = $a18;
	}
	elsif ($cov == 19)
	{
	    $a{$chrpos} = $a19;
	}
	elsif ($cov == 20)
	{
	    $a{$chrpos} = $a20;
	}
	elsif ($cov == 21)
	{
	    $a{$chrpos} = $a21;
	}
	elsif ($cov == 22)
	{
	    $a{$chrpos} = $a22;
	}
	elsif ($cov == 23)
	{
	    $a{$chrpos} = $a23;
	}
	elsif ($cov == 24)
	{
	    $a{$chrpos} = $a24;
	}
	else
	{
	    print "Wrong\n";
	}
    }
}
close FH;

%pi = ();
%theta = ();
%nan = ();

open FH, "<SN_M_SNP.txt.sites.pi";
while (<FH>)
{
    if (/(\S+)\s+(\S+)\s+(\S+)/)
    {
	$chr = $1;
	$pos = $2;
	$p = $3;
	$chrpos = $chr . "_" . $pos;
	if ($p eq "-nan")
	{
	  $nan{$chrpos} = 0;
	}
	elsif ($p > 0)
	{
	    $pi{$chrpos} = $p;
	    $theta{$chrpos} = 1 / $a{$chrpos};
	}
    }
}
close FH;
close OUT;

open FH, "<SN_M_10TO150_count_tab_025cov_pos.depth";
open OUT, ">SN_M_10TO150_count_tab_025cov_pos_diversity.depth";
$n = 0;
print OUT "Num\tChr\tPos\ttheta\tpi\n";
while (<FH>)
{
    $n = $n + 1;
    if (/(\S+)\s+(\S+)/)
    {
	$chr = $1;
	$pos = $2;
	$chrpos = $chr . "_" . $pos;
	if (!(exists ($nan{$chrpos})))
	{
	    if (exists ($pi{$chrpos}))
	    {
		print OUT "$n\t$chr\t$pos\t$theta{$chrpos}\t$pi{$chrpos}\n";
	    }
	    else
	    {
		print OUT "$n\t$chr\t$pos\t0\t0\n";
	    }
	}
    }
}
close FH;
close OUT;
