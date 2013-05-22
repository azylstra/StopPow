package stoppowgui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;
import javax.swing.JMenuItem;
import javax.swing.JPopupMenu;

/**
 * Implement a run method for updating the GUI. This class is used
 * to launch updates to the panel GUI through use of the Event Dispatch Thread,
 * for example via:
 * {@code SwingUtilities.invokeLater(new CalculatorPanelGUIRunnable(this));}
 * @class CalculatorPanelGUIRunnable
 * @author Alex Zylstra
 * @date 2013/05/18
 */
class CalculatorPanelGUIRunnable implements Runnable
{
    /** The CalculatorPanel in question */
    private CalculatorPanel window;
    
    /** 
     * Create a new Runnable object for the given panel.
     * @param window the Panel for GUI updates
     */
    public CalculatorPanelGUIRunnable( CalculatorPanel window )
    {
        this.window = window;
    }

    /**
     * Calls the CalculatorPanel updateGUI method. Designed to be used
     * with SwingUtilities.invokeLater to launch via Event Dispatch Thread
     */
    @Override
    public void run() {
        window.updateGUI();
    }
    
}

/**
 * Implement generic functionality for a JPanel taking input from the
 * user and displaying a result for stopping power calculations.
 * This class defines common methods for all CalculatorPanels, and
 * also adds functionality for adjusting output formatting and
 * dealing with parts of the treading infrastructure.
 * @class CalculatorPanel
 * @author Alex Zylstra
 * @date 2013/05/18
 */
abstract class CalculatorPanel extends javax.swing.JPanel implements Runnable
{
    /** The key for the model to be used in the calculation. */
    protected String key;
    /** The ModelManager in use by this application */
    protected ModelManager manager;
    /** A ThreadPool for calculations */
    protected ThreadPoolExecutor exec;
    
    /** Number of digits to display */
    protected int displayedDigits;
    /** Default number of digits to display */
    protected final static int DEFAULT_DISPLAYED_DIGITS = 2;
    /** A String description of the format for displaying output */
    protected String outputFormat;
    /** Default (ie initial) format string */
    protected final static String defaultFormat = "%.2f";
    /** Which mode to use for output (decimal or exponential) */
    protected int outputMode;
    /** display as floating point */
    protected final static int OUTPUT_MODE_FP = 0; 
    /** display in exponential notation */
    protected final static int OUTPUT_MODE_EXP = 1;
    
    /** Allow for locking of the top-level panel infrastructure */
    protected final Lock panelLock = new ReentrantLock();
    
    /** A popup menu which appears over calculated results to allow
     * for changing the formatting.
     */
    protected JPopupMenu rightClickMenu;
    
    /**
     * Create a new CalculatorPanel.
     * @param key The String key for the model that this panel is using. 
     * Must match the key for model in the ModelManager. Can be null
     * when instantiating a new panel.
     * @param manager The ModelManager in use by this application
     * @param calcPool A ThreadPool for use when doing background calculations
     */
    public CalculatorPanel(String key, ModelManager manager, ThreadPoolExecutor calcPool)
    {
        // set some class members:
        this.key = key;
        this.manager = manager;
        exec = calcPool;
        
        // set default result formatting:
        outputFormat = defaultFormat;
        displayedDigits = DEFAULT_DISPLAYED_DIGITS;
        outputMode = OUTPUT_MODE_FP;
        generateFormatString();
        
        // set up the right click menu:
        // this is to allow the user to change some formatting stuff
        // for displayed outputs
        rightClickMenu = new JPopupMenu();
        
        // item to set decimal format:
        JMenuItem decimal = new JMenuItem("Decimal");
        decimal.addActionListener(new ActionListener(){
            @Override
            public void actionPerformed(ActionEvent evt) {
                displayDecimal();
                updateGUI();
            }
            
        });
        rightClickMenu.add(decimal);
        
        // item to set exponential format:
        JMenuItem exponential = new JMenuItem("Scientific");
        exponential.addActionListener(new ActionListener(){
            @Override
            public void actionPerformed(ActionEvent evt) {
                displayExponential();
                updateGUI();
            }
            
        });
        rightClickMenu.add(exponential);
        
        // item to add displayed digits
        JMenuItem more = new JMenuItem("More");
        more.addActionListener(new ActionListener(){
            @Override
            public void actionPerformed(ActionEvent evt) {
                displayMoreDigits();
                updateGUI();
            }
            
        });
        rightClickMenu.add(more);
        
        // item to remove displayed digits
        JMenuItem less = new JMenuItem("Less");
        less.addActionListener(new ActionListener(){
            @Override
            public void actionPerformed(ActionEvent evt) {
                displayFewerDigits();
                updateGUI();
            }
            
        });
        rightClickMenu.add(less);
    }
    
    /**
     * Update the underlying stopping power calculation. Does not launch in
     * a new thread, use the Runnable interface for that:
     * {@code exec.execute(myCalculatorPanelShift)}
     * where exec is a ThreadPool.
     */
    abstract public void updateCalculation();
    
    /**
     * Clear the displayed result.
     */
    abstract public void clearResults();
    
    /**
     * Update the displayed dynamic elements to reflect the latest
     * calculation result and mode selection.
     */
    abstract public void updateGUI();
      
    /**
     * Update the stopping power model being used for calculations
     * in this window.
     * @param newKey the String key used in the ModelManager for newModel.
     */
    public void updateModel(String newKey)
    {
        panelLock.lock();
        try
        {
            key = newKey;
        }
        finally
        {
            panelLock.unlock();
        }
        
        // launch calculation update in new thread:
        exec.execute(this);
    }
    
    /**
     * Run the calculation. For use with threads.
     */
    @Override
    public void run() 
    {
        updateCalculation();
    }
    
    // ---------------------------------------
    //  Adjust display of calculated results
    // ---------------------------------------
    /**
     * Display one additional digit in the output.
     */
    public void displayMoreDigits()
    {
        displayedDigits++;
        generateFormatString();
    }
    
    /**
     * Display one fewer digit in the output (minimum of 1).
     */
    public void displayFewerDigits()
    {
        // make sure we always have at least 1 displayed digit:
        if( displayedDigits > 1 )
        {
            displayedDigits--;
            generateFormatString();
        }
    }
    
    /**
     * Toggle display to be decimal.
     */
    public void displayDecimal()
    {
        outputMode = OUTPUT_MODE_FP;
        generateFormatString();
    }
    
    /**
     * Toggle display to be exponential (scientific).
     */
    public void displayExponential()
    {
        outputMode = OUTPUT_MODE_EXP;
        generateFormatString();
    }
    
    /**
     * Generate the format string for the current parameters.
     */
    private void generateFormatString()
    {
        String newFormat = "%.";
        
        newFormat += Integer.toString(displayedDigits);
        
        if( outputMode == OUTPUT_MODE_EXP )
            newFormat += "e";
        else
            newFormat += "f";
        
        outputFormat = newFormat;
    }
}
