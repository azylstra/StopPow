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

import cStopPow.FloatVector2D;
import cStopPow.StopPow;
import cStopPow.cStopPow;
import javax.swing.SwingUtilities;

/**
 * Generate a Eout vs Ein/Thickness given Thickness/Ein plot.
 * @author Alex Zylstra
 */
public class PlotGenerator_Eout extends PlotGenerator {
    
    /**
     * Create a new PlotGenerator object.
     * @param models The ModelManager in use by this application.
     * @param param the parameters to use when generating the plot.
     */
    public PlotGenerator_Eout(ModelManager models, PlotParameters param)
    {
        super(models,param);
    }

    /**
     * Get the ordinate (x axis) label for this plot.
     * @return The label to use as a String
     */
    @Override
    protected String xLabel() 
    {
        String xLabel = "";
        if( param.secondary.contains("Ein") )
            xLabel = "Thickness (" + getModeUnits() + ")";
        else
            xLabel = "Ein (MeV)";
        return xLabel;
    }
    
    /**
     * Get the abscissa (y axis) label for this plot.
     * @return The label to use as a String
     */
    @Override
    protected String yLabel()
    {
        return "Eout (MeV)";
    }


    /**
     * Create a single dataset for Eout plots
     * @param key the key corresponding to the StopPow model to use
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
        // secondary is either "Ein" or "Thickness"
        if( param.secondary.contains("Ein") )
        {
            cStopPow.get_Eout_vs_Thickness(model,param.secondaryValue,data);
        }
        else
        {
            cStopPow.get_Eout_vs_Ein(model,param.secondaryValue,data);
        }

        return createFloatArray(data);
    }
}
