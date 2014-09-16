// StopPow - a charged-particle stopping power library
// Copyright (C) 2014  Massachusetts Institute of Technology / Alex Zylstra

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

package stoppowgui;

import SciTK.DialogError;
import SciTK.PlotXYLine;
import cStopPow.FloatVector2D;
import cStopPow.StopPow;
import cStopPow.cStopPow;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import javax.swing.SwingUtilities;

/**
 * Wrap parameters for plot generation into a class.
 * @author Alex Zylstra
 * @date 2013/06/05
 */
class PlotParameters
{
    public String[] keys;
    public String abscissa;
    public String ordinate;
    public String secondary;
    public float secondaryValue;
    public String mode;
}


/**
 * Generate plots from defined parameters.
 * This class implements Runnable so that it can be
 * invoked in a Thread, e.g.
 * {@code exec.execute(gen);}
 * where exec is a ThreadPoolExecutor and gen is a PlotGenerator.
 * The behavior when invoked like that (or via run()) is
 * to create the plot and display it.
 * @author Alex Zylstra
 * @date 2013/06/05
 */
abstract class PlotGenerator implements Runnable
{
    protected ModelManager manager;
    protected PlotParameters param;


    /**
     * Create a new PlotGenerator object.
     * @param models The ModelManager in use by this application.
     * @param param the parameters to use when generating the plot.
     */
    public PlotGenerator(ModelManager models, PlotParameters param)
    {
        this.manager = models;
        this.param = param;
    }
    
    /**
     * Generate and display all plots.
     */
    @Override
    public void run() 
    {
        // go through plots and get a list of all keys used
        ArrayList<String> modelKeys = new ArrayList<String>();
        modelKeys.addAll(Arrays.asList(param.keys));
                
        // first, get the model locks:
        for(String key : modelKeys)
            manager.lockModel(key);

        // wrap rest of the code in try/finally
        try
        {
            createPlot();
        }
        finally
        {
            // release all locks:
            for(String key : modelKeys)
                manager.unlockModel(key);
        }
    }
    
    protected abstract float[][] createDataset(String key);
    protected abstract String xLabel();
    protected abstract String yLabel();
    
    /**
     * Create a single plot from parameters
     */
    private void createPlot()
    {
        // single plot only:
        if( param.keys.length == 1 )
        {
            // get the data to plot:
            float[][] dataset = createDataset(param.keys[0]);
            
            // and show it:
            SwingUtilities.invokeLater(new PlotXYLineGenerator(dataset,param.keys[0],xLabel(),yLabel()));
        }
        // multi plot:
        else if( param.keys.length > 1 )
        {
            // data for the plot:
            float[][][] dataset = new float[param.keys.length][1][1];
            
            // construct all data sets:
            for(int i=0; i<param.keys.length; i++)
            {
                dataset[i] = createDataset(param.keys[i]);
            }
            
            // make the plot:
            SwingUtilities.invokeLater(new PlotXYLineGenerator(dataset,param.keys,xLabel(),yLabel()));            
        }
    }
    
    protected float[][] createFloatArray(FloatVector2D input)
    {
        float[][] ret;
        ret = new float[(int)input.size()][(int)input.get(0).size()];
        
        // iterate over input to populate ret:
        for(int i=0; i<input.size(); i++)
        {
            for(int j=0; j<input.get(i).size(); j++)
                ret[i][j] = input.get(i).get(j);
        }
        
        return ret;
    }
    
    
    /**
     * Get a string describing the units for length quantities
     * @return "um" or "mg/cm2" depending on the mode selection
     */
    public String getModeUnits()
    {
        // set the unit strings
        String units = param.mode;
        
        // remove MeV/ and parens:
        units = units.replace("MeV/","");
        units = units.replace("(","");
        units = units.replace(")","");
        units = units.trim();
        
        return units;
    }
}

/**
 * Generate a PlotXYLine. This class is written to allow
 * invocation via the Event Dispatch Thread when the calculation
 * is running in another thread, e.g. via
 * {@code SwingUtilities.invokeLater(new PlotXYLineGenerator(datasets,keys));}
 * @author Alex Zylstra
 * @date 2013/05/21
 */
class PlotXYLineGenerator implements Runnable
{
    /** The datasets to use for the plot */
    private float[][][] data;
    /** The names of above datasets */
    private String[] name;
    /** An ordinate label for the plot */
    private String xLabel;
    /** An abscissa label for the plot */
    private String yLabel;

    /**
     * Construct a PlotXYLine from a single set of data.
     * Note: Plot is not generated until run() is invoked.
     * @param data a 2xn array of values to plot
     * @param name The name of the dataset (e.g. for legend)
     * @param xLabel a label for the ordinate (x axis)
     * @param yLabel a label for the abscissa (y axis)
     */
    public PlotXYLineGenerator(float[][] data, String name, String xLabel, String yLabel)
    {
        this.data = new float[1][data.length][data[0].length];
        this.data[0] = data;
        this.name = new String[1];
        this.name[0] = name;
        
        this.xLabel = xLabel;
        this.yLabel = yLabel;
    }

    /**
     * Construct a PlotXYLine from multiple sets of data.
     * Note: Plot is not generated until run() is invoked.
     * @param data a mx2xn array (m datasets, n points each)
     * @param name The names of the datasets (e.g. for legend)
     * @param xLabel a label for the ordinate (x axis)
     * @param yLabel a label for the abscissa (y axis)
     */
    public PlotXYLineGenerator(float[][][] data, String[] name, String xLabel, String yLabel)
    {
        this.data = data;
        this.name = name;
        
        this.xLabel = xLabel;
        this.yLabel = yLabel;
    }

    /**
     * Invoke this method to generate and display the plot.
     */
    @Override
    public void run() {
        try
        {
            PlotXYLine p = new PlotXYLine(data,name,xLabel,yLabel);
            p.show();
        }
        catch(Exception e)
        {
            DialogError msg = new DialogError(null,e.getMessage());
        }
    }
}