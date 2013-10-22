<? 
	$Title="Figures"; include "../header.php"; 

	// Tour listing
	function begin_tours($title, $link)
	{ 
		echo '<div class="post"><div class="header"><h3><a name="' . $link . '"></a>' . $title . '</h3></div><div class="content"><p><ol>';
	}
	function end_tours()
	{  
		echo '</ol></p></div></div>';
	}
	function tour($link, $title)
	{
		echo '<li><a href="' . $link . '/">' . $title . '</a></li>';
	}
	
	// Table of content
	function begin_toc()
	{ 
	 	echo '<p align="center"><div align="center"><font size="+1"><b>'; 
	 	echo '- Table of contents -'; 
	 	echo '</b></font></div>'; 
	 	echo '<p align="center"><div align="center"><font size="+0">';
	}
	function end_toc()
	{
	 	echo '</font></div><br/><br/>';
	}
	function toc_entry($title, $link)
	{ 
		echo '<a href="#' . $link . '">' . $title . '</a><br/>';
	}
	
	// What's news section
	function begin_news()
	{
		echo '<div class="post"><div class="header"><h3> What\'s new</h3></div><div class="content"><p><ul>';
	}
	function end_news()
	{
		echo '</ul></p></div></div>';
	}
	function news_display($date, $content, $link)
	{
		echo '<li>(' . $date . ') ' . $content . ' [<a href="' . $link . '">link</a>]' . '</li>';
	}
?>
	
	

<?
	include "index_news.php";
	include "index_tours.php";
	include "../footer.php"; 
?>