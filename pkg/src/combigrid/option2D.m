X = [ 0.9 , 0.935086 , 0.96366 , 0.981203 , 0.990267 , 0.994611 , 0.996882 , 0.998496 , 1 , 1.0015 , 1.00312 , 1.00539 , 1.00973 , 1.0188 , 1.03634 , 1.06491 , 1.1]; 
 Y = [ 0.9 , 0.96366 , 0.990267 , 0.996882 , 1 , 1.00312 , 1.00973 , 1.03634 , 1.1]; 
 res = [ 0.000854479 , 0.00393591 , 0.00794085 , 0.00934287 , 0.0100321 , 0.0107229 , 0.0122275 , 0.0181266 , 0.0433539 ; 0.00183042 , 0.00795873 , 0.0148494 , 0.0168773 , 0.0178297 , 0.0187639 , 0.0206772 , 0.0285616 , 0.0639188 ; 0.00435218 , 0.014601 , 0.0235442 , 0.0259004 , 0.0269154 , 0.0278789 , 0.0297908 , 0.0393467 , 0.0804927 ; 0.00633053 , 0.0209354 , 0.0313197 , 0.0340637 , 0.0353331 , 0.0365825 , 0.0391828 , 0.0502712 , 0.0973789 ; 0.00758347 , 0.0245319 , 0.0357344 , 0.0386676 , 0.0400444 , 0.0414154 , 0.0443032 , 0.0567767 , 0.106458 ; 0.0081739 , 0.0261648 , 0.0378808 , 0.0409152 , 0.0423388 , 0.0437685 , 0.046793 , 0.0597445 , 0.110826 ; 0.00848228 , 0.0269963 , 0.0390028 , 0.0421002 , 0.0435556 , 0.0450139 , 0.0481005 , 0.06137 , 0.113027 ; 0.00868233 , 0.0275719 , 0.0398013 , 0.0429443 , 0.0444226 , 0.0458998 , 0.0490274 , 0.0625093 , 0.114578 ; 0.0088681 , 0.0281034 , 0.0405447 , 0.0437316 , 0.045231 , 0.0467269 , 0.0498928 , 0.0635385 , 0.116061 ; 0.00905529 , 0.0286288 , 0.0412872 , 0.0445194 , 0.0460403 , 0.0475548 , 0.0507593 , 0.0645548 , 0.117548 ; 0.00925584 , 0.0291764 , 0.0420817 , 0.0453648 , 0.0469103 , 0.0484455 , 0.0516876 , 0.0656419 , 0.119148 ; 0.00954136 , 0.0299339 , 0.0431969 , 0.0465551 , 0.0481396 , 0.0497102 , 0.052995 , 0.0671665 , 0.121396 ; 0.0100648 , 0.0313337 , 0.0453344 , 0.0488621 , 0.0505211 , 0.0521547 , 0.0555516 , 0.070071 , 0.125793 ; 0.011216 , 0.0340141 , 0.0500003 , 0.0538498 , 0.0556521 , 0.0574191 , 0.0611026 , 0.0764832 , 0.134395 ; 0.014534 , 0.0406746 , 0.059874 , 0.0641546 , 0.0661523 , 0.0681163 , 0.0722544 , 0.0901635 , 0.150339 ; 0.0226613 , 0.0547815 , 0.0786281 , 0.0831093 , 0.0853088 , 0.0874897 , 0.0923631 , 0.118603 , 0.181483 ; 0.0281138 , 0.0797295 , 0.0932777 , 0.101255 , 0.104919 , 0.108584 , 0.116325 , 0.149802 , 0.213628]; 
 [x,y]=meshgrid(Y,X);
 surf(x,y,res); 
 