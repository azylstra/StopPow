/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package stoppowgui;

import SciTK.*;

import javax.swing.JToggleButton;
import javax.swing.JLabel;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import javax.swing.AbstractButton;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.SwingUtilities;


/**
 *
 * @author alex
 */
public class StopPowGUI extends SciTK_App
{
    private ModelManagerGUI models;
    
    /** Constructor */
    public StopPowGUI()
    {
        // use constructor for superclass SciTK_App with specified size:
        super(200,200);
        NAME = "Stopping Power Calculator";
        AUTHOR = "Alex Zylstra";
        DATE = "May 07, 2013";
        
        //String logo_path = "resources/logo.png";
        //about_icon = new ImageIcon(getClass().getResource(logo_path));
        about_icon = new ImageIcon();
        //tray_icon = new ImageIcon(getClass().getResource(logo_path));
        tray_icon = new ImageIcon();
        
        //initUI();
        this.setSize(200, 200);
        setLayout(new BoxLayout(getContentPane(), BoxLayout.Y_AXIS));
	setDefaultCloseOperation(EXIT_ON_CLOSE);
        app_ui();
        
        // set up model manager:
        models = new ModelManagerGUI(this);
    }
    
    /** Handle opening a file. Does nothing for this application. */
    public void open_file() {}
    
    /** Load the initial UI */
    public void app_ui()
    {
        // Set a bit of text:
        JLabel instruction = new JLabel(" Stopping power: ");
        instruction.setFont(new Font("Serif", Font.BOLD, 16));
        instruction.setAlignmentX(0.5f);
        add(instruction);

        // model_manager_ button
        JToggleButton button_model_manager = new JToggleButton("Configure Models");
        button_model_manager.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                // find out if the button is selected:
                AbstractButton button = (AbstractButton) event.getSource();
                boolean selected = button.getModel().isSelected();
                
                if( selected )
                {
                    // selected, show the model manager:
                    models.setVisible(true);
                    button.setSelected(true);
                    
                } 
                else if( !selected )
                {
                    // deselected, hide it:
                    models.setVisible(false);
                    button.setSelected(false);
                }
            }
        });
        button_model_manager.setAlignmentX(0.5f);
        add(button_model_manager);
        
        // Calculators button
        JToggleButton button_calculators = new JToggleButton("Calculators");
        button_calculators.addActionListener(new ActionListener() {
            public void itemStateChanged(ItemEvent event) {
                // foo
            }

            @Override
            public void actionPerformed(ActionEvent ae) {
                // don't need to do anything
            }
        });
        button_calculators.setAlignmentX(0.5f);
        add(button_calculators);
        
        // plots button
        JToggleButton button_plots = new JToggleButton("Plots");
        button_plots.addActionListener(new ActionListener() {
            public void itemStateChanged(ItemEvent event) {
                // foo
            }

            @Override
            public void actionPerformed(ActionEvent ae) {
                // don't need to do anything
            }
        });
        button_plots.setAlignmentX(0.5f);
        add(button_plots);


    }
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // First try to load the native library directly from a file:
        try
        {
            System.loadLibrary("cStopPow"); 
        }
        catch( UnsatisfiedLinkError e )
        {
            //SciTK_Dialog d = new Dialog_Error("Could not load native library.");
            //return;
        }
        
        // the native library might also be packaged:
        try
        {
            //StopPowLibLoader s = new StopPowLibLoader();
            //s.loadStopPow();
            NativeUtils.loadLibraryFromJar("/cStopPow/libcStopPow.jnilib");
        }
        catch( UnsatisfiedLinkError e )
        {
            SciTK_Dialog d = new Dialog_Error("Could not load native library.");
            return;
        }
        catch( java.io.IOException e )
        {
            SciTK_Dialog d = new Dialog_Error(" IO Error: " + e.getMessage());
        }
        
        // Create a new NouveauAnalysis object and let
        // Java invoke it when appropriate:
        SwingUtilities.invokeLater(new Runnable() {
                public void run() {
                        StopPowGUI g = new StopPowGUI();
                        g.setVisible(true);
                }
        });
    }
}
