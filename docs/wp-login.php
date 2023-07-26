<!DOCTYPE html>
	<html xmlns="http://www.w3.org/1999/xhtml" lang="en-US">
	<head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title>rDock &rsaquo; Log In</title>
	<link rel='stylesheet' id='wp-admin-css'  href='http://rdock.sourceforge.net/wp-admin/css/wp-admin.min.css?ver=3.7.1' type='text/css' media='all' />
<link rel='stylesheet' id='buttons-css'  href='http://rdock.sourceforge.net/wp-includes/css/buttons.min.css?ver=3.7.1' type='text/css' media='all' />
<link rel='stylesheet' id='colors-fresh-css'  href='http://rdock.sourceforge.net/wp-admin/css/colors-fresh.min.css?ver=3.7.1' type='text/css' media='all' />
<meta name='robots' content='noindex,nofollow' />
	</head>
	<body class="login login-action-login wp-core-ui">
	<div id="login">
		<h1><a href="http://wordpress.org/" title="Powered by WordPress">rDock</a></h1>
	
<form name="loginform" id="loginform" action="http://rdock.sourceforge.net/wp-login.php" method="post">
	<p>
		<label for="user_login">Username<br />
		<input type="text" name="log" id="user_login" class="input" value="" size="20" /></label>
	</p>
	<p>
		<label for="user_pass">Password<br />
		<input type="password" name="pwd" id="user_pass" class="input" value="" size="20" /></label>
	</p>
		<p class="forgetmenot"><label for="rememberme"><input name="rememberme" type="checkbox" id="rememberme" value="forever"  /> Remember Me</label></p>
	<p class="submit">
		<input type="submit" name="wp-submit" id="wp-submit" class="button button-primary button-large" value="Log In" />
		<input type="hidden" name="redirect_to" value="http://rdock.sourceforge.net/wp-admin/" />
		<input type="hidden" name="testcookie" value="1" />
	</p>
</form>

<p id="nav">
	<a href="http://rdock.sourceforge.net/wp-login.php?action=lostpassword" title="Password Lost and Found">Lost your password?</a>
</p>

<script type="text/javascript">
function wp_attempt_focus(){
setTimeout( function(){ try{
d = document.getElementById('user_login');
d.focus();
d.select();
} catch(e){}
}, 200);
}

wp_attempt_focus();
if(typeof wpOnload=='function')wpOnload();
</script>

	<p id="backtoblog"><a href="http://rdock.sourceforge.net/" title="Are you lost?">&larr; Back to rDock</a></p>
	
	</div>

	
		<div class="clear"></div>
	</body>
	</html>
	