
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="http://<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://r-forge.r-project.org/"><img src="http://<?php echo $themeroot; ?>/imagesrf/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->

<p> The main focus of the R Package Fitob is the PDE-based financial derivative pricing.  
Fitob is able to solve not just the Black-Scholes PDE in a multidimensional setting, 
but also the more general Fokker-Planck (convection-diffusion) PDE. 
Fitob provides a general scripting interface that can describe almost any financial 
contract in an efficient and unique way. <br> 
The following examples illustrate the PDE-based pricing approach of financial contracts:
</p>
<p> <img src="heston.png" alt="Pricing a 1D Option under the Heston process" width="304" height="228"> </p>

<p> <img src="sim_S_A.gif" alt="Early withdrawal boundary for a Guaranteed Minimum Withdrawal Benefit (GMWB). 
A is the state of the account, S is the underlying's price" width="304" height="228"> </p>

<p> Manual and Tutorial of the Fitob R Package <strong> <a href="http://fitob.r-forge.r-project.org/fitob_manual.pdf">  here </strong> </a> </p>

<p> Download Fitob R Package <strong> <a href="http://r-forge.r-project.org/R/?group_id=1515">  here </strong> </a> </p>

<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

</body>
</html>
