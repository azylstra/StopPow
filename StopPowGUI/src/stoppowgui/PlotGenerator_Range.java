package stoppowgui;

import cStopPow.FloatVector2D;
import cStopPow.StopPow;
import cStopPow.cStopPow;

/**
 * Generate a range vs energy in plot.
 * @author Alex Zylstra
 */
public class PlotGenerator_Range extends PlotGenerator {
    
    /**
     * Create a new PlotGenerator object.
     * @param models The ModelManager in use by this application.
     * @param param the parameters to use when generating the plot.
     */
    public PlotGenerator_Range(ModelManager models, PlotParameters param)
    {
        super(models,param);
    }

    /**
     * Create a single dataset for range plots
     * @param key the key corresponding to the StopPow model to use
     * @param mode The calculation mode ("MeV/um" or "MeV/(mg/cm2)") to use
     * @return The data as 2xn float array
     */
    @Override
    protected float[][] createDataset(String key)
    {
        // get the model:
        StopPow model = manager.get_model(key);

        // return value:
        FloatVector2D data = new FloatVector2D();

        // set mode appropriately:
        if( param.mode.equals("MeV/um") )
            model.set_mode( StopPow.getMODE_LENGTH() );
        else if( param.mode.equals("MeV/(mg/cm2)") )
            model.set_mode( StopPow.getMODE_RHOR() );

        // build array:
        cStopPow.get_Range_vs_E(model, data);

        return createFloatArray(data);
    }

    /**
     * Get the ordinate (x axis) label for this plot.
     * @return The label to use as a String
     */
    @Override
    protected String xLabel()
    {
        return "Energy (MeV)";
    }

    /**
     * Get the abscissa (y axis) label for this plot.
     * @return The label to use as a String
     */
    @Override
    protected String yLabel()
    {
        return "Thickness (" + getModeUnits() + ")";
    }
}
