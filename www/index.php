
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
<!--
<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?> -->

<!-- end of project description -->

<center> <h1>Welcome to the Fitob R Package<br></h1> </center>

<p> The main focus of the R Package Fitob is the PDE-based financial derivative pricing.  
Fitob is able to solve not just the Black-Scholes PDE in a multi-dimensional setting, 
but also the more general Fokker-Planck (convection-diffusion) PDE. 
Fitob provides a general scripting interface that can describe almost any financial 
contract in an efficient and unique way. <br> </p>
<br>
<br>

<p>
Important links:
<ul>
<li>Manual and Tutorial of the Fitob R Package <strong> <a href="http://fitob.r-forge.r-project.org/fitob_manual.pdf">  here </strong> </a> </li>
<li>Download Fitob R Package <strong> <a href="http://r-forge.r-project.org/R/?group_id=1515">  here </strong> </a> </li>
<li>The project summary page you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a> </li>
</ul>
</p>
<br>
<br>

<p>
The following examples illustrate the PDE-based pricing approach of financial contracts (you can find more examples 
with detailed explanations in the examples section of the 
<strong> <a href="http://fitob.r-forge.r-project.org/fitob_manual.pdf"> manual </strong> </a>):
</p>
<br>
<br>

   <ul>
   <h5>
   Pricing a 1D American Option. The plot shows the resulting price curve. <br>
   <center>
   <img src="american1D.png" alt="Pricing a 1D Option under the Heston process" width="504" height="428" align="middle"> <br>
   </center>
   You can download the necessary files under:
   </strong> <a href="http://fitob.r-forge.r-project.org/examples/american1D.xml">  american1D.xml </strong> </a> ,
   </strong> <a href="http://fitob.r-forge.r-project.org/examples/american1DScript">  american1DScript </strong> </a> <br>
   In order to run the example, install Fitob and type in R:<br>
   >  require(fitob); <br>
   >  liMy = fitobPriceMesh("american1D.xml" , "american1DScript",5);  <br>
   >  fitobMeshPlot(liMy); <br>
   >  liMy[[1]] #the scalar price<br>
   </h5>

   <ul>
   <h5>
   Pricing a 1D Option under the Heston process. The plot shows the resulting price surface in 2D. <br>
   <center>
   <img src="heston.png" alt="Pricing a 1D Option under the Heston process" width="504" height="428" align="middle"> <br>
   </center>
   You can download the necessary files under:
   </strong> <a href="http://fitob.r-forge.r-project.org/examples/Heston2D.xml">  Heston2D.xml </strong> </a> , 
   </strong> <a href="http://fitob.r-forge.r-project.org/examples/Heston2D_Script">  Heston2D_Script </strong> </a>
   In order to run the example, install Fitob and type in R:<br>
   >  require(fitob); <br>
   >  liMy = fitobPriceMesh("Heston2D.xml" , "Heston2D_Script",5);  <br>
   >  fitobMeshPlot(liMy); <br>
   >  liMy[[1]] #the scalar price <br>
   </h5>

   </ul>
   <h5>
   Early withdrawal boundary for a Guaranteed Minimum Withdrawal Benefit (GMWB) modelled in 5D. 
   A is the account of the benefit, S is the underlying's price. (see J.Benk, D. Plueger: Hybrid Parallel Solutions of the Black-Scholes PDE with the Truncated Combination Technique. In Proceedings of the HPCS conference, 2012 Madrid.) <br>
   <center>
   <img src="sim_S_A.gif" alt="Early withdrawal boundary for a Guaranteed Minimum Withdrawal Benefit (GMWB). 
   A is the state of the account, S is the underlying's price" width="504" height="428" align="middle"> <br>
   </center>
   </h5>

</body>
</html>
