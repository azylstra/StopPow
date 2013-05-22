package stoppowgui;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
 
import SciTK.DialogError;
import java.util.ArrayList;

/**
 * Loader for native libraries. Automatically tries the current classpath,
 * and then if necessary loads from a JAR file.
 * 
 * 
 * @brief Load a native library
 * @class LibraryLoader
 * @author Alex Zylstra
 * @date 2013/05/10
 */
public class LibraryLoader {
 
    /** Constructor. All methods are static, so not really needed */
    private LibraryLoader() {}
 
    /**
     * Load a native library.
     * @param name the name of the library you want to load
     * @return whether the library was loaded successfully
     */
    public static boolean loadLibrary(String name)
    {        
        // first, try a normal load procedure:
        try
        {
            System.loadLibrary(name); 
        }
        catch(UnsatisfiedLinkError e)
        {
            return false;
        }
        
        // next, try to load from JAR:
        try
        {
            load_JAR(name);
        }
        catch(IllegalArgumentException e)
        {
            DialogError err = new DialogError(null,"Couldn't open library: " + e.getMessage());
            return false;
        }
        catch(IOException e)
        {
            DialogError err = new DialogError(null,"Couldn't open library: " + e.getMessage());
            return false;
        }
        catch(UnsatisfiedLinkError e)
        {
            DialogError err = new DialogError(null,"Couldn't open library: " + e.getMessage());
            return false;
        }
        
        return true;
    }
    
    /**
     * Loads library from current JAR archive
     * 
     * The file from JAR is copied into system temporary directory and then loaded. The temporary file is deleted after exiting.
     * Method uses String as filename because the pathname is "abstract", not system-dependent.
     * 
     * @param filename The filename inside JAR as absolute path (beginning with '/'), e.g. /package/File.ext
     * @throws IOException
     * @throws IllegalArgumentException
     */
    public static void load_JAR(String name) throws IOException 
    {
        // sanity check:
        if (!name.startsWith("/"))
            name = "/" + name;
        
        // parse the name:
        String[] parts = parse_name(name);
 
        // Prepare temporary file
        File temp = File.createTempFile(parts[0], parts[2]);
        temp.deleteOnExit();
 
        if (!temp.exists())
            throw new FileNotFoundException("File " + temp.getAbsolutePath() + " does not exist.");
 
        // Prepare buffer for data copying
        byte[] buffer = new byte[1024];
        int readBytes;
 
        // Open and check input stream
        InputStream is = LibraryLoader.class.getResourceAsStream(name);
        if (is == null) {
            throw new FileNotFoundException("File " + name + " was not found inside JAR.");
        }
 
        // Open output stream and copy data between source file in JAR and the temporary file
        OutputStream os = new FileOutputStream(temp);
        try {
            while ((readBytes = is.read(buffer)) != -1) {
                os.write(buffer, 0, readBytes);
            }
        } finally {
            // If read/write fails, close streams safely before throwing an exception
            os.close();
            is.close();
        }
 
        // Finally, load the library
        try
        {
            System.load(temp.getAbsolutePath());
        }
        catch(Exception e)
        {
            DialogError err = new DialogError(null,"Couldn't open library: " + e.getMessage());
        }
    }
    
    /**
     * Split a filename into a prefix (path to file), name of the file itself,
     * and suffix (i.e. extension).
     * @param name the String to split
     * @return parts a String[] of 3 elements containing prefix, name, suffix
     */
    private static String[] parse_name(String name)
    {
        // Obtain filename from path
        String[] parts = name.split("/");
        String filename = null;
        // if there are multiple parts:
        if( parts.length > 1 )
            filename = parts[parts.length-1];
 
        // Split filename to prexif and suffix (extension)
        String prefix = "";
        String suffix = null;
        if (filename != null) 
        {
            parts = filename.split("\\.", 2);
            prefix = parts[0];
            if( parts.length > 1)
                suffix = "." + parts[parts.length-1];
        }
 
        // Check if the filename is okay
        if (filename == null || prefix.length() < 3)
        {
            throw new IllegalArgumentException("The filename has to be at least 3 characters long.");
        }
        
        // build return array:
        String[] ret = new String[3];
        ret[0] = prefix;
        ret[1] = filename;
        ret[2] = suffix;
        
        return ret;
    }
}
