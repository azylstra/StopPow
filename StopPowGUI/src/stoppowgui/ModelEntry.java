package stoppowgui;

import cStopPow.StopPow;

/**
 * Wrapper class for model data stored in the ModelManager
 * @class ModelEntry
 * @date 2013/05/10
 * @author Alex Zylstra
 */
public class ModelEntry extends java.lang.Object
{
    /**
     * Constructor
     * @param name_in name (key) of the model
     * @param type_in type of model (e.g. SRIM)
     * @param model_in the StopPow object
     * @param info_in a general info string (can be anything)
     */
    public ModelEntry(String name_in, String type_in, StopPow model_in, String info_in)
    {
        super();
        name = name_in;
        model = model_in;
        type = type_in;
        info = info_in;
    }
    
    /** Name (key) of the model */
    public String name;
    /** Type of model (e.g. SRIM) */
    public String type;
    /** The StopPow object */
    public StopPow model;
    /** A general info String (can be anything) */
    public String info;
}
