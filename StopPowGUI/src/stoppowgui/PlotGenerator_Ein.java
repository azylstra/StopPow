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

/**
 * Generate a Ein vs Thickness/Eout given Eout/Thickness plot.
 * @author Alex Zylstra
 */
public class PlotGenerator_Ein extends PlotGenerator {
    
    /**
     * Create a new PlotGenerator object.
     * @param models The ModelManager in use by this application.
     * @param param the parameters to use when generating the plot.
     */
    public PlotGenerator_Ein(ModelManager models, PlotParameters param)
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
        String xLabel;
        if( param.secondary.contains("Eout") )
            xLabel = "Thickness (" + getModeUnits() + ")";
        else
            xLabel = "Eout (MeV)";
        return xLabel;
    }
    
    /**
     * Get the abscissa (y axis) label for this plot.
     * @return The label to use as a String
     */
    @Override
    protected String yLabel()
    {
        return "Ein (MeV)";
    }

    /**
     * Create a single dataset for Ein plots
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

        // secondary is either "Eout" or "Thickness"
        if( param.secondary.contains("Eout") )
        {
            cStopPow.get_Ein_vs_Thickness(model, param.secondaryValue, data);
        }
        else
        {
            cStopPow.get_Ein_vs_Eout(model, param.secondaryValue, data);
        }

        return createFloatArray(data);
    }
}
