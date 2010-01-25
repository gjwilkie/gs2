require 'lib/autodoc.rb'

# Get sourceforge username
username = FileTest.exist?('username.txt') ? File.read('username.txt').chomp : (puts "Please enter your sourceforge username (you must have admin privileges for gyrokinetics. To avoid seeing this message please run 'echo <your username> > username.txt'"; STDIN.gets.chomp)

######################################
# Setup autodoccer
autodoccer = Autodoc.new(Dir.pwd + '/../trunk', Dir.pwd + '/gs2_documentation', {:sub_directories => ['utils', 'geo']})
autodoccer.code_name = "GS2"
autodoccer.code_website = "http://gyrokinetics.sourceforge.net"
autodoccer.code_description = "GS2 is a gyrokinetic flux tube initial value turbulence code which can be used for fusion or astrophysical plasmas."
autodoccer.custom_tabs = {"Wiki" => "http://sourceforge.net/apps/mediawiki/gyrokinetics", "Old Website" => "http://gs2.sourceforge.net/"}
autodoccer.ignore_files = ['test_os', 'utils/redistribute.f90', 'gs2', 'ingen', 'rungridgen', 'fortdep']
autodoccer.welcome_message = %[<p>Welcome to the GS2 Documentation Home Page. This is the documentation that comes from comments made in the source code. For other help see the gyrokinetics <a href="http://sourceforge.net/apps/mediawiki/gyrokinetics">wiki</a> or the old <a href="http://gs2.sourceforge.net/">website</a>.</p>

<p>To add to this documentation, just surround useful comments in the source code by &lt;doc&gt;&lt;/doc&gt; braces. Comments within a subroutine or function definition document that subroutine or function, and comments between the beginning of a module and the first subroutine document the module. Anything that is surrounded by double square brackets: [[name]] will link to the wiki page of the same name. See <a href="file_utils.html">file_utils</a> as an example of a documented module.</p> 

<p>To update this html (if you are impatient to see your new comments up here), install <a href = "http://www.andre-simon.de/zip/download.html">Highlight</a> then run &quot;ruby document_gs2.rb&quot  in the folder gyrokinetics/gs2/autodoc (requires <a href = "http://www.ruby-lang.org/en/downloads/">Ruby</a> version 1.9 or higher).<p>]
######################################


######################################
# Optional - used mainly for debugging Autodoc - can be commented out
autodoccer.strict = true
autodoccer.external_modules = ['hdf5']
autodoccer.external_globals = ['fftw_f77_create_plan', 'rfftwnd_f77_create_plan', 'rfftw2d_f77_create_plan', 'rfftw_f77_destroy_plan', 'h5pcreate_f', 'h5pset_dxpl_mpio_f', 'h5pset_fapl_mpio_f', 'h5pclose_f', 'h5screate_simple_f', 'h5sselect_hyperslab_f','h5sclose_f', 'h5dwrite_f', 'h5dcreate_f', 'h5screate_f', 'h5dclose_f', 'h5open_f', 'h5dopen_f', 'h5fcreate_f', 'h5fclose_f', 'h5close_f', 'pgctab', 'pgline', 'pgqls', 'pgqci', 'pgsls', 'pgsci', 'pgsvp', 'pgswin', 'pgbox', 'vtsymdef', 'vtbegin', 'vtend', 'mpi_initialized', 'mpi_init', 'mpi_comm_size', 'mpi_comm_rank', 'mpi_finalize', 'mpi_bcast', 'mpi_reduce', 'mpi_allreduce', 'mpi_barrier', 'mpi_send', 'mpi_ssend', 'mpi_recv', 'mpi_cart_create', 'mpi_cart_coords', 'mpi_cart_sub']
#######################################

autodoccer.produce_highlighted_source_code = true
autodoccer.document_fortran_intrinsic = false #true #Documents fortran instrinsic functions. Makes autodoc run more slowly.

#######################################
# This section links any [[name]] surrounded by double square brackets to the
# wiki article with the same name
autodoccer.customize_documentation = Proc.new do |highlighted|
	highlighted.gsub(/\[\[(\w+)\]\]/, '<a href="http://sourceforge.net/apps/mediawiki/gyrokinetics/index.php?title=\1">\1</a>')
end
autodoccer.customize_highlighting = Proc.new do |highlighted|
	highlighted.gsub(/\[\[((?:\s*\<\/span\>\s*)?)(\w+)((?:\s*\<span[^>]*\>\s*)?)\]\]/, '[[\1<a href="http://sourceforge.net/apps/mediawiki/gyrokinetics/index.php?title=\2">\2</a>\3]]')
end
#######################################

 
autodoccer
autodoccer.write_documentation

string = "rsync -av --delete gs2_documentation/  #{username},gyrokinetics@web.sourceforge.net:htdocs/gs2_documentation/"

puts string
exec string



# $not_found.uniq!
# Dir.chdir('../temp') do
# 	$not_found.each do |function|
# # 		puts function
# 		lines = `grep #{function} */* *`.split("\n")
# # 		puts lines
# 		puts lines.grep(Regexp.new("\\s#{Regexp.escape(function)}\\s")).grep(/^[^!]*(subroutine|interface)/)
# 	end
# end
# p $not_found.sort
