package stoppowgui;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ThreadPoolExecutor;

/**
 *
 * @author alex
 */
public class PlotManager {
    private ThreadPoolExecutor exec;
    private PlotManagerGUI parent;
    private ModelManager models;
    private Map<Integer,PlotGenerator> plots;
    private Integer count;
    private List<Integer> queue;
    
    public PlotManager(ThreadPoolExecutor exec, ModelManager models)
    {
        this.exec = exec;
        this.models = models;
        plots = new HashMap<Integer,PlotGenerator>() {};
        queue = new ArrayList<Integer>() {};
        count = 0;
    }
    
    /**
     * Add a plot to the generator. However, the plot is not generated
     * until this object's run method is invoked.
     * @param newParam A PlotParameters object describing the plot to be created.
     */
    public void addPlot(PlotParameters newParam)
    {
        PlotGenerator g = null;
        // add based on type:
        if( newParam.abscissa.equals("dE/dx") )
            g = new PlotGenerator_dEdx(models,newParam);
        else if( newParam.abscissa.equals("Ein") )
            g = new PlotGenerator_Ein(models,newParam);
        else if( newParam.abscissa.equals("Eout") )
            g = new PlotGenerator_Eout(models,newParam);
        else if( newParam.abscissa.equals("Thickness") )
            g = new PlotGenerator_Thickness(models,newParam);
        else if( newParam.abscissa.equals("Range") )
            g = new PlotGenerator_Range(models,newParam);
        else
            return; // base case, did not find type, return and do not add
        
        plots.put(count,g);
        queue.add(count);
        count++;
    }
    
    /**
     * Generate the plots
     */
    public void generate()
    {
        PlotGenerator newPlot;
        for(Integer i : queue)
        {
            newPlot = plots.get(i);
            exec.execute(newPlot);
        }
        queue.clear();
    }
}
