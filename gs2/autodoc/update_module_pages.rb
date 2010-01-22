#!/home/edmundhighcock/Tools/bin/bin/ruby

require 'cgi'
require 'fileutils'
require 'pp'


class Autodoc
	class ParsingError < StandardError
	end
	attr_accessor :code_website, :code_name, :code_description, :welcome_message, :produce_highlighted_source_code, :external_modules, :external_globals, :strict
	attr_reader :function_references, :modules
	def initialize(source_dir, html_dir, options={})
		@source_dir, @html_dir = source_dir, html_dir
		@function_references = {} #eval(File.read("function_references.rb"))
		@documentation_references = {}
		@out_pages = {}
		@highlighted_files = {}
		@sub_directories = options[:sub_directories] ? options[:sub_directories] + ['.'] : ['.']
		@title_prefix = (options[:title_prefix] or "Autodoc Documentation")
		@modules = {}
		@uses = {}
		@external_modules = []
		@external_globals = []
		@strict = false
	end
	def exclude?(file)
		return (not File.file? file or ['Makefile', 'README', 'test_os'].include? file or file =~ /Makefile/ or ['.inc', '.svn', '.in'].include? File.extname(file) or file =~ /\~$/ or (false and FileTest.exist?(out_file) and  File.mtime(file) < File.mtime(out_file)))
	end
# 	def analyse_and_highlight_source_code
# 		analyse_documentation
# 		Dir.chdir(@source_dir) do
# 		@sub_directories.each do |dir| 
# 			Dir.chdir(dir) do
# 				Dir.entries(Dir.pwd).each do |file|	
# 					sub_dir = dir == "." ? "" : dir + "/"		
# 					out_file = "#{sub_dir}#{file}.html"
# 					(next) if exclude? file
# #					puts file
# 					@highlighted_files[@html_dir + "/source/" + out_file] =  analyse_and_highlight_file(file, out_file)
# 				end
# 			end
# 
# 		end
# 
# 
# 	end
# end

# File.open("function_references.rb", 'w'){|ref_file| ref_file.puts function_references.pretty_inspect} 

#puts function_references.pretty_inspect

#function_references.dup.each do |function, reference|
#	puts reference
#	function_references[function] = reference.gsub(/(?<double>[^\/]+\/)/){puts $~; gets; $~[:double]}
#end

#http://gyrokinetics.sourceforge.net/autodoc/

# 	end #analyse_and_highlight_source_code
	FORTRAN_INTRINSIC = ["I", "ABORT", "ABS", "ACCESS", "ACHAR", "ACOS", "ACOSH", "ADJUSTL", "ADJUSTR", "AIMAG", "AINT", "ALARM", "ALL", "ALLOCATED", "AND", "ANINT", "ANY", "ASIN", "ASINH", "ASSOCIATED", "ATAN", "ATAN", "ATANH", "BESSEL", "BESSEL", "BESSEL", "BESSEL", "BESSEL", "BESSEL", "BIT", "BTEST", "C", "C", "C", "C", "C", "C", "CEILING", "CHAR", "CHDIR", "CHMOD", "CMPLX", "COMMAND", "COMPLEX", "CONJG", "COS", "COSH", "COUNT", "CPU", "CSHIFT", "CTIME", "DATE", "DBLE", "DCMPLX", "DFLOAT", "DIGITS", "DIM", "DOT", "DPROD", "DREAL", "DTIME", "EOSHIFT", "EPSILON", "ERF", "ERFC", "ERFC", "ETIME", "EXIT", "EXP", "EXPONENT", "FDATE", "FLOAT", "FGET", "FGETC", "FLOOR", "FLUSH", "FNUM", "FPUT", "FPUTC", "FRACTION", "FREE", "FSEEK", "FSTAT", "FTELL", "GAMMA", "GERROR", "GETARG", "GET", "GET", "GETCWD", "GETENV", "GET", "GETGID", "GETLOG", "GETPID", "GETUID", "GMTIME", "HOSTNM", "HUGE", "HYPOT", "IACHAR", "IAND", "IARGC", "IBCLR", "IBITS", "IBSET", "ICHAR", "IDATE", "IEOR", "IERRNO", "INDEX", "INT", "INT", "INT", "IOR", "IRAND", "IS", "IS", "ISATTY", "ISHFT", "ISHFTC", "ISNAN", "ITIME", "KILL", "KIND", "LBOUND", "LEADZ", "LEN", "LEN", "LGE", "LGT", "LINK", "LLE", "LLT", "LNBLNK", "LOC", "LOG", "LOG", "LOG", "LOGICAL", "LONG", "LSHIFT", "LSTAT", "LTIME", "MALLOC", "MATMUL", "MAX", "MAXEXPONENT", "MAXLOC", "MAXVAL", "MCLOCK", "MCLOCK", "MERGE", "MIN", "MINEXPONENT", "MINLOC", "MINVAL", "MOD", "MODULO", "MOVE", "MVBITS", "NEAREST", "NEW", "NINT", "NOT", "NULL", "OR", "PACK", "PERROR", "PRECISION", "PRESENT", "PRODUCT", "RADIX", "RAN", "RAND", "RANDOM", "RANDOM", "RANGE", "REAL", "RENAME", "REPEAT", "RESHAPE", "RRSPACING", "RSHIFT", "SCALE", "SCAN", "SECNDS", "SECOND", "SELECTED", "SELECTED", "SELECTED", "SET", "SHAPE", "SIGN", "SIGNAL", "SIN", "SINH", "SIZE", "SIZEOF", "SLEEP", "SNGL", "SPACING", "SPREAD", "SQRT", "SRAND", "STAT", "SUM", "SYMLNK", "SYSTEM", "SYSTEM", "TAN", "TANH", "TIME", "TIME", "TINY", "TRAILZ", "TRANSFER", "TRANSPOSE", "TRIM", "TTYNAM", "UBOUND", "UMASK", "UNLINK", "UNPACK", "VERIFY", "XOR"]
	def highlight_source
		each_source_file do |file, subdir_with_slash|
			@highlighted_files["#@html_dir/source/#{subdir_with_slash}#{file}.html"] =  highlight_file(file)
		end
		exit
	end
	def highlight_file(file) 
		case File.extname(file)
		when "fpp"
			syntax = " -s=f90"
		else
			syntax = " "
		end
# 		return %x[highlight -H -a #{syntax} -i #{file}  --style lucretia --inline-css -K 12 -k Monaco -l].gsub(/(\d+\s*\<\/span\>\s*(?:\<span[^>]*\>)?\s*subroutine\s*(?:\<\/span\>\s*)?(?:\<span[^>]*\>\s*)?)(\w+)(\s*.{40})/m){$2; @function_references[$2] = "#{out_file}\##$2"; %[#$1<a name="#$2"></a>#$2#$3]}
		highlighted =  %x[highlight -H -a #{syntax} -i #{file}  --style greenlcd --inline-css -K 12 -k Monaco -l]
		module_blocks = highlighted.scan(/.+?(?:\d+\s*\<\/span\>\s*(?:\<span[^>]*\>)?\s*module(?:\s*\<\/span\>)?(?![^\n]*procedure)|\Z)/m)
# 		p module_blocks.size
		highlighted = module_blocks.shift
		module_blocks.each do |block|
			(highlighted += block; next) if block.length == 1 #The last character
# 			p block.length
# 			p block.split("\n")[0]
# 			p block[0, 8]
# 			p block.scan(/\A\s*(\S+)/)
# 			p 'hello'
# 			
# 			
			module_name = block.scan(/\A\s*(\S+)/)[0][0]
# 			p module_name

			subroutine_blocks = block.scan(/.+?(?:\d+\s*\<\/span\>\s*(?:\<span[^>]*\>)?\s*subroutine(?:\s*\<\/span\>)?(?:\s*\<span[^>]*\>)?|\Z)/m)
			block = subroutine_blocks.shift
			subroutine_blocks.each do |subblock|
				(block += subblock; next) if subblock.length == 1 #The last character
				subroutine_name = subblock.scan(/\A\s*(\w+)/)[0][0]
				(p subblock; exit) if subroutine_name =~ /span/
				subblock.gsub!(/^[^!]*(call(?:(?:(?:\s*\<\/span\>\s*)?(?:\s*\<span[^>]*\>\s*)?)|(?:\s+)))(\w+)/m) do
					before = $1; call_name = $2
					next if @external_globals.include? call_name 
					used = nil
					@uses[[module_name, subroutine_name]].each do |use_mod, use_subroutines|
						if use_subroutines.include? call_name
							used = [use_mod, call_name]
						end
					end if @uses[[module_name, subroutine_name]]
					next if used and  @external_modules.include? used[0] 
					reference = (@function_references[used] or @function_references[[module_name, call_name]])
					reference = "http://gcc.gnu.org/onlinedocs/gcc-4.4.3/gfortran/#{call_name.upcase}.html##{call_name.upcase}" if not reference and FORTRAN_INTRINSIC.include? call_name.upcase
# 					'reference', reference
					(p 'used', used, 'file', file, 'module_name', module_name, 'subname', subroutine_name, 'call_name', call_name; raise ParsingError) if @strict and not reference
					reference ? %[#{before}<a href= "#{reference}">#{call_name}</a>] : "#{before}#{call_name}"
				end
				block += subblock
			end
			block.gsub!(/(\d+\s*\<\/span\>\s*(?:\<span[^>]*\>)?\s*subroutine\s*(?:\<\/span\>\s*)?(?:\<span[^>]*\>\s*)?)(\w+)(\s*.{40})/m){%[#$1<a name="#$2"></a>#$2#$3]}
# 			puts block
# 			puts "\n\n\n"
			highlighted += block
		end
# 		p @function_references
# 		puts highlighted
# 		exit
# 		
# 		.gsub(/(call(?:(?:(?:\s*\<\/span\>\s*)?(?:\s*\<span[^>]*\>\s*)?)|(?:\s+)))(\w+)/m){@function_references[$2] ? %[#$1<a href= "#{@function_references[$2]}">#$2</a>] : "#$1#$2"} 
	end
	def write_documentation
		analyse_source unless @analysed_source
		highlight_source if @produce_highlighted_source_code
		@highlighted_files.each do |out_file, text|
			puts out_file
			File.open(out_file, 'w'){|write_file| write_file.puts text}
		end
		Dir.chdir(@html_dir) do 
			@modules.each do |module_name, data|
				File.open("#{module_name}.html", 'w')do |file| 
					mod = ModulePage.new(module_name, data, self)
					file.puts mod
				end
			end
			File.open("index.html", 'w')do |file| 
				index = IndexPage.new(self)
				file.puts index
			end
			File.open("moduleindex.html", 'w')do |file| 
				index = ModuleIndexPage.new(self)
				file.puts index
			end


		end
	end #def write_documentation

	def each_source_file(&block)
		Dir.chdir(@source_dir) do
		@sub_directories.each do |dir| 
			Dir.chdir(dir) do
				Dir.entries(Dir.pwd).each do |file|	
					subdir_with_slash = dir == "." ? "" : dir + "/"		
# 					out_file = "#{sub_dir}#{file}.html"
					(next) if exclude? file
					yield(file, subdir_with_slash)
				end # Dir.entreis
			end #Dir.chdir
		end # @sub_directories.each
		end # Dir.chdir(source_dir)		
	end
	
	#not File.file? file or ['Makefile', 'README', 'test_os'].include? file or file =~ /Makefile/ or ['.inc', '.svn', '.in'].include? File.extname(file) or file =~ /\~$/ or (false and FileTest.exist?(out_file) and  File.mtime(file) < File.mtime(out_file)) 
	
	def analyse_source
		@analysed_source = true
		each_source_file do |file, subdir_with_slash|
			begin
				text = File.read(file)
# 				next unless text  =~ /\<wkdoc\>/
			rescue
				next
			end
			modules = text.scan(/(^\s*module\s*(?!procedure)(\S+).*?\s*end\s+module)/m)
# 					pp modules[0]
			modules.each do |modtext, name|
				(p file, name, modtext; raise Autodoc::ParsingError) if name =~ /procedure/
				analyse_module(name, file, subdir_with_slash, modtext)
			end
		end
		File.open('uses.rb', 'w'){|file| file.puts @uses.pretty_inspect}
		File.open('fr.rb', 'w'){|file| file.puts @function_references.pretty_inspect}

	end # analyse_source
	def analyse_module(modname, file, subdir_with_slash, modtext)
		@modules[modname] = {}
		@modules[modname][:file_name] = file
		@modules[modname][:subroutines] = {}
		if @produce_highlighted_source_code
			interfaces = modtext.scan(/^\s*interface\s*(\w+)/)
			(p modname, file, interfaces; raise ParsingError) if (["interface", "contains"] - interfaces).size < 2
			interfaces.flatten.each do |interface|
				@function_references[[modname, interface]] = %[#{subdir_with_slash}#{file}.html##{modname}_#{interface}] 
			end
		end
		subroutine_blocks = modtext.split(/^\s*subroutine\s*/)
		@modules[modname][:description] = subroutine_blocks[0].scan(/\<wkdoc\>(.*?)\<\/wkdoc\>/m).map{|match| match[0].gsub(/\n\s*\!/, '')}.join("\n")
		subroutine_blocks.slice(1..(subroutine_blocks.size)).each do |block|
			name = block.scan(/(\A\w+)/)[0][0]
# 						p block
# 			puts "\n\n\n"
			function_call = block.scan(/(\A.+(?:\&\s)?.+)/)[0][0].sub(/[\&\n]/, '')
# 						function_call = block.scan(/(\A.*?(?:\(.*[\n]?.*\))?)/)[0][0].sub(/[\&\n]/, '')
			comments = block.scan(/\<wkdoc\>(.*?)\<\/wkdoc\>/m).map{|match| match[0]}
			comments = comments.map do |comment|
				comment.gsub(/\n\s*\!/, '')
			end
# 			p 
			use_statements = block.scan(/(^\s*use[^&\n\r]+(?:(?:\&\s+)?[^&\n\r]+)+)/)
# 			(p use_statements[0][0]; exit) if use_statements.find{|statements| statements[0] =~ /\&/} 
			@uses[[modname, name]] = {}
			use_statements.each do |statement_list|
				words = statement_list[0].gsub(/\&/, '').scan(/\w+/)
				words.delete("use"); words.delete("only")
				use_mod = words.shift
				use_subroutines = words
				@uses[[modname, name]][use_mod] ||= []
				@uses[[modname, name]][use_mod] += words
				@uses[[modname, name]][use_mod].uniq!
# 				p words; exit
			end
			@modules[modname][:subroutines][name] = [function_call, comments]
			@function_references[[modname, name]] = %[#{subdir_with_slash}#{file}.html##{modname}_#{name}] if @produce_highlighted_source_code
			@documentation_references[[modname, name]] = %[#{modname}.html##{name}]
		end
	end
	class Page
		def to_s
			<<EOF
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="content-type" content="text/html; charset=utf-8" />
<title>Metamorphosis Design Free Css Templates</title>
<meta name="keywords" content="" />
<meta name="description" content="" />
<link href="styles.css" rel="stylesheet" type="text/css" />
</head>
<body>
<div id="main"><!-- start header -->
<div id="header">
 <div id="logo">
	<h1>#{top_label}</h1>
	<h2><a href="http://www.metamorphozis.com/" id="metamorph">Designed by Metamorphosis Design,</a></h2><h2><a href="http://autodoc.rubyforge.org/" id="metamorph">produced by Autodoc</a></h2>
	</div>
	<div id="menu">
		<ul>
			<li><a href="index.html">Home</a></li>
			<li><a href="moduleindex.html">Modules</a></li>
			<!--<li><a href="https://sourceforge.net/apps/mediawiki/gyrokinetics">Wiki</a></li>
			<li><a href="#">About</a></li>
			<li><a href="#">Contact</a></li>-->
		</ul>
	</div>	
<!-- end header -->
</div>
<!-- start page -->
<div id="top"></div>	
<div id="page">
#{title_box}
#{content}
	<!-- start sidebar -->
	<div id="sidebar2" class="sidebar">
			#{subroutine_sidebar}
			#{module_sidebar}
	</div>
	<!-- end sidebar -->
<div style="clear: both;">&nbsp;</div>

</div>
	<div id="bottom"></div>
<!-- end page -->
<!-- start footer -->
<div id="footer">
	 <p>Copyright &copy; 2009. <a href="#">Privacy Policy</a> | <a href="#">Terms of Use</a> | <a href="http://validator.w3.org/check/referer" title="This page validates as XHTML 1.0 Transitional"><abbr title="eXtensible HyperText Markup Language">XHTML</abbr></a> | <a href="http://jigsaw.w3.org/css-validator/check/referer" title="This page validates as CSS"><abbr title="Cascading Style Sheets">CSS</abbr></a></p> 
	<p>Design by <a href="http://www.metamorphozis.com/" title="Free Flash Templates">Free Flash Templates</a>
		</p>
</div>
</div>
<!-- end footer -->
</body>
</html>

EOF
		end # def page
		def top_label
			%[<a href = "#{@autodoccer.code_website or "#"}">#{@autodoccer.code_name} Documentation</a>]
		end
		def module_sidebar
				keys = @autodoccer.modules.keys
			<<EOF
				<h2>Module Index</h2>
				<ul class="back_title">
				#{keys.slice(1..keys.size).inject(%[\n\t\t\t<li class="top"><a href="#{keys[0]}.html">#{keys[0]}</a></li>]) do |str, name|
					str + %[\n\t\t\t<li><a href="#{name}.html">#{name}</a></li>]
				end}
				</ul>
EOF
		end #def module_sidebar
		def subroutine_sidebar
			return ""
		end #def subroutine_sidebar

	end # class page
	class ModulePage < Page
		def initialize(module_name, data, autodoccer)
			@module_name, @module_file_name, @subroutines, @autodoccer = module_name, data[:file_name], data[:subroutines], autodoccer
			@description = data[:description]
			 @function_references = autodoccer.function_references
		end
		def subroutine_div(name, function_call, comments)
		lines = comments.map do |comment|
			line = "<li>#{comment}</li>"
			@function_references.each do |(modname, name), reference|
				next unless modname == name
				line.gsub(/name/, 
					%[<a href="#{reference}">name</a>])
			end
			line
		end
		return <<EOF
			#{@function_references[[modname, name]] ? %[<h2 class="title"><a href="source/#{@function_references[[modname, name]]}" name="#{name}">#{name}</a></h2>] : %[<h2 class="title">#{name}</h2>]} 
			
		<div class="entry"><small>Call Prototype:</small> #{function_call} #{@function_references[[modname, name]] ? %[<small><a href="source/#{@function_references[[modname, name]]}">View Source</a></small>] : %[]} </div>
			<ul>
				#{lines.join("\n\t\t\t")}
			</ul><br>
EOF
		
		end #subroutine_div
		def content
		<<EOF
	<!-- start content -->
	<div id="content">
#{@subroutines.inject("") do |str, (name, (function_call, comments))|
			str + subroutine_div(name, function_call, comments)
end}
	</div>
	<!-- end content -->
EOF
		end # def content
		def title_box
			<<EOF
<div id="box">
	<h1><a href="#">Module: #@module_name</a></h1>
<p>#{@description or "Brief description coming soon!"}</p><br>
<div class="entry"><small>Last updated #{Time.now.to_s}</small></div>
</div>
EOF
		end #def title_box
		def subroutine_sidebar
				keys = @subroutines.keys
			<<EOF
				<h2>Subroutines</h2>
				<ul class="back_title">
				#{keys.slice(1..keys.size).inject(%[\n\t\t\t<li class="top"><a href="##{keys[0]}">#{keys[0]}</a></li>]) do |str, name|
					str + %[\n\t\t\t<li><a href="##{name}">#{name}</a></li>]
				end}
				</ul>
EOF
		end #def subroutine_sidebar


	end # class ModulePage
	
	class IndexPage < Page
		def initialize(autodoccer)
			@autodoccer =autodoccer
			@function_references = autodoccer.function_references
		end
		def welcome_message
			return <<EOF
			<h2 class="title">Welcome!</h2>
			<p>#{@autodoccer.welcome_message}</p>
			<div class="entry"><small> This documentation has been automatically generated from the #{@autodoccer.code_name} source code by Autodoc</small></div>
EOF
		end # def welcome_message
		def content
		<<EOF
	<!-- start content -->
	<div id="content">
#{welcome_message}
	</div>
	<!-- end content -->
EOF
		end # def content
		def title_box
			<<EOF
<div id="box">
	<h1><a href="#">#{@autodoccer.code_name} Documentation Home Page</a></h1>
<p>#{@autodoccer.code_description}</p><br>
<div class="entry"><small>Last updated #{Time.now.to_s}</small></div>
</div>
EOF
		end #def title_box


	end # class IndexPage

	class ModuleIndexPage < Page
		def initialize(autodoccer)
			@autodoccer =  autodoccer
			@modules = {}
			@autodoccer.modules.each do |name, data|
				@modules[name] = [data[:file_name], data[:description]]
			end
			 @function_references = autodoccer.function_references
		end
		def module_div(name, file_name, description)
# 			p file_name
		return <<EOF
<h2 class="title"><a href="#{name}.html" name="#{name}">#{name}</a></h2>
		<div class="entry"><small>Source File:</small><a href = "source/#{file_name}.html"> #{file_name}</a></div>
		<div class="entry">#{description}</div><br>
EOF
		
		end #module_div
		def content
		<<EOF
	<!-- start content -->
	<div id="content">
#{@modules.inject("") do |str, (name, (file_name, description))|
			str + module_div(name, file_name, description)
end}
	</div>
	<!-- end content -->
EOF
		end # def content
		def title_box
			<<EOF
<div id="box">
	<h1><a href="#">Modules</a></h1>
<p>A list of modules with short descriptions.</p><br>
<div class="entry"><small>Last updated #{Time.now.to_s} using Autodoc</small></div>
</div>
EOF
		end #def title_box
		def subroutine_sidebar
			return ""
		end #def subroutine_sidebar


	end # class ModuleIndexPage

end #class Autodoc

autodoccer = Autodoc.new(Dir.pwd + '/../temp', Dir.pwd + '/documentation', {:sub_directories => ['utils', 'geo']})
autodoccer.code_name = "GS2"
autodoccer.code_website = "http://gyrokinetics.sourceforge.net"
autodoccer.code_description = "GS2 is a gyrokinetic flux tube initial value turbulence code which can be used for fusion or astrophysical plasmas."
autodoccer.strict = true
autodoccer.external_modules = ['hdf5']
autodoccer.external_globals = ['fftw_f77_create_plan', 'rfftwnd_f77_create_plan', 'rfftw2d_f77_create_plan', 'rfftw_f77_destroy_plan']

# autodoccer.analyse_and_highlight_source_code
autodoccer.produce_highlighted_source_code = true
autodoccer.write_documentation