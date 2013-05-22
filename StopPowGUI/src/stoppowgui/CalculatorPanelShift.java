package stoppowgui;

import cStopPow.StopPow;
import java.awt.Color;
import java.awt.Font;
import java.util.concurrent.ThreadPoolExecutor;
import javax.swing.SwingUtilities;


/**
 * Generate a JPanel containing UI elements allowing the user
 * to do calculations for energy shift. The user can select
 * calculations for energy downshift, upshift, or material
 * thickness for a given dE.
 * @author Alex Zylstra
 * @class CalculatorPanelShift
 * @date 2013/05/18
 */
public class CalculatorPanelShift extends CalculatorPanel
{
    // for display of text:
    /** The default Font for display of UI text labels. */
    private Font defaultFont;
    /** The Font to use when displaying special text,
     * ie the selected calculation mode. */
    private Font selectedFont;
    /** The default Color for displaying UI text labels. */
    private Color defaultColor;
    /** The Color to use when displaying special text,
     * ie the selected calculation mode. */
    private Color selectedColor;
    
    // for keeping track of the calculation mode:
    /** The current calculation mode */
    private int calcMode;
    /** Static member to denote upshift calculations */
    private final static int MODE_EIN = 0;
    /** Static member to denote thickness calculations */
    private final static int MODE_THICKNESS = 1;
    /** Static member to denote downshift calculations */
    private final static int MODE_EOUT = 2;
    
    /** The calculation result */
    private float calcResult;
    
    /**
     * Create a new CalculatorPanelShift.
     * @param key The String key for the model that this panel is using. 
     * Must match the key for model in the ModelManager. Can be null
     * when instantiating a new panel.
     * @param manager The ModelManager in use by this application
     * @param calcPool A ThreadPool for use when doing background calculations
     */
    public CalculatorPanelShift(String key, ModelManager manager, ThreadPoolExecutor calcPool)
    {
        super(key,manager,calcPool);
        
        initComponents();
        
        // set fonts (piggyback on GUI):
        defaultFont = EinLabel.getFont();
        selectedFont = EoutLabel.getFont();
        defaultColor = EinLabel.getForeground();
        selectedColor = EoutLabel.getForeground();
        
        // start in Eout mode:
        calcMode = MODE_EOUT;
        
        // start w/o a result:
        calcResult = Float.NaN;
    }
    
    /**
     * Update the underlying stopping power calculation. Does not launch in
     * a new thread, use the Runnable interface for that:
     * {@code exec.execute(myCalculatorPanelShift)}
     * where exec is a ThreadPool.
     */
    @Override
    public void updateCalculation() 
    {
        // sanity check:
        if( key == null || !manager.containsKey(key) )
        {
            calcResult = Float.NaN;
            return;
        }
        
        // do a quasi lock on the input fields:
        try {
            EinTextField.setEditable(false);
            ThicknessTextField.setEditable(false);
            EoutTextField.setEditable(false);
        }
        catch( Exception e ) {}
        
        // try to lock the panel and the model:
        panelLock.lock();
        manager.lockModel(key);
        try
        {
            // get the model
            StopPow model = manager.get_model(key);
                    
            // get the values currently displayed in three fields
            // N.B. only 2 will be used
            float Ein = Float.NEGATIVE_INFINITY, Thickness=Float.NEGATIVE_INFINITY, Eout=Float.NEGATIVE_INFINITY;
            // only parse if text strings are not empty:
            // wrap in try/catch in case text is empty or otherwise non-numeric
            try {
                Ein = Float.parseFloat( EinTextField.getText() );
            }
            catch(java.lang.NumberFormatException e) {}
            try {
                Thickness = Float.parseFloat( ThicknessTextField.getText() );
            }
            catch(java.lang.NumberFormatException e) {}
            try {
                Eout = Float.parseFloat( EoutTextField.getText() );
            }
            catch(java.lang.NumberFormatException e) {}

            // set the model's mode based on drop-down menu:
            int modeIndex = modeComboBox.getSelectedIndex();
            if( modeComboBox.getItemAt(modeIndex) == "MeV/(mg/cm2)")
                model.set_mode( StopPow.getMODE_RHOR() );
            else if( modeComboBox.getItemAt(modeIndex) == "MeV/um")
                model.set_mode( StopPow.getMODE_LENGTH() );
                        
            // do appropriate calculation for current mode,
            // and set the displayed text for it
            // wrap in try/catch for exceptions thrown by native code
            try {
                switch(calcMode)
                {
                    case MODE_EIN:
                        calcResult = model.Ein(Eout, Thickness);
                        break;
                    case MODE_THICKNESS:
                        calcResult = model.Thickness(Ein, Eout);
                        break;
                    // Calculate energy out from energy in and thickness
                    case MODE_EOUT:
                        calcResult = model.Eout(Ein, Thickness);
                        break;
                }
            }
            catch (java.lang.IllegalArgumentException e) {
                calcResult = Float.NaN;
            }
        }
        finally
        {
            panelLock.unlock();
            manager.unlockModel(key);
        }
        
        // need to update the GUI elements
        // use SwingUtilities to let the EDT do this later:
        SwingUtilities.invokeLater(new CalculatorPanelGUIRunnable(this));
        
    }
    
    /**
     * Clear the displayed result.
     */
    @Override
    public void clearResults()
    {
        switch(calcMode)
        {
            case MODE_EIN: 
                EinTextField.setText("");
                break;
            case MODE_THICKNESS:
                ThicknessTextField.setText("");
                break;
            case MODE_EOUT:
                EoutTextField.setText("");
                break;
        }
    }
    
    /**
     * Choose the mode to be used from user input.
     * For CalculatorPanelShift, this is done by clicking on the
     * text labels, therefore this takes a {@link MouseEvent}.
     * @param evt 
     */
    private void chooseMode(java.awt.event.MouseEvent evt)
    {
        // set the mode based on the source of the event:
        if( evt.getSource() == EinLabel )
        {            
            calcMode = MODE_EIN;
        }
        else if( evt.getSource() == ThicknessLabel )
        {            
            calcMode = MODE_THICKNESS;
        }
        else if( evt.getSource() == EoutLabel )
        {            
            calcMode = MODE_EOUT;
        }
        
        // immediately ask Swing to update the GUI:
        SwingUtilities.invokeLater(new CalculatorPanelGUIRunnable(this));
        // update calculation:
        exec.execute(this);
        // the above will also update GUI after calculation is done
    }
    
    /**
     * Update the displayed dynamic elements to reflect the latest
     * calculation result and mode selection.
     */
    @Override
    public void updateGUI()
    {
        // first set defaults for all things:
        // default is to let text fields be edited:
        EinTextField.setEditable(true);
        ThicknessTextField.setEditable(true);
        EoutTextField.setEditable(true); 
        // set default fonts:
        EinLabel.setFont(defaultFont);
        ThicknessLabel.setFont(defaultFont);
        EoutLabel.setFont(defaultFont);
        // and label text colors:                   
        EinLabel.setForeground(defaultColor);            
        ThicknessLabel.setForeground(defaultColor);
        EoutLabel.setForeground(defaultColor);
        
        // Now set the appropriate value and special properties
        // depending on the mode:
        switch(calcMode)
        {
            case MODE_EIN:          
                if( Float.isNaN(calcResult) )
                    EinTextField.setText( "" );
                else
                    EinTextField.setText( String.format(outputFormat, calcResult) );
                EinTextField.setEditable(false);
                EinLabel.setFont(selectedFont);
                EinLabel.setForeground(selectedColor);                                    
                break;
            case MODE_THICKNESS:    
                if( Float.isNaN(calcResult) )
                    ThicknessTextField.setText( "" );
                else
                    ThicknessTextField.setText( String.format(outputFormat, calcResult) );
                ThicknessTextField.setEditable(false);
                ThicknessLabel.setFont(selectedFont);
                ThicknessLabel.setForeground(selectedColor);                                    
                break;
            case MODE_EOUT:     
                if( Float.isNaN(calcResult) )
                    EoutTextField.setText( "" );
                else
                    EoutTextField.setText( String.format(outputFormat, calcResult) );    
                EoutTextField.setEditable(false);     
                EoutLabel.setFont(selectedFont);
                EoutLabel.setForeground(selectedColor);                               
                break;
        }
    }
    
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        EinLabel = new javax.swing.JLabel();
        EinTextField = new javax.swing.JTextField();
        jLabel2 = new javax.swing.JLabel();
        jLabel3 = new javax.swing.JLabel();
        modeComboBox = new javax.swing.JComboBox();
        EoutLabel = new javax.swing.JLabel();
        EoutTextField = new javax.swing.JTextField();
        ThicknessLabel = new javax.swing.JLabel();
        ThicknessTextField = new javax.swing.JTextField();
        jLabel6 = new javax.swing.JLabel();
        jLabel7 = new javax.swing.JLabel();

        EinLabel.setText("E in:");
        EinLabel.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                EinLabelMouseClicked(evt);
            }
        });

        EinTextField.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseReleased(java.awt.event.MouseEvent evt) {
                EinTextFieldMouseReleased(evt);
            }
        });
        EinTextField.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                EinTextFieldActionPerformed(evt);
            }
        });

        jLabel2.setText("MeV");

        jLabel3.setText("Mode: ");
        jLabel3.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                jLabel3MouseClicked(evt);
            }
        });

        modeComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "MeV/um", "MeV/(mg/cm2)" }));
        modeComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                modeComboBoxActionPerformed(evt);
            }
        });

        EoutLabel.setFont(new java.awt.Font("Cantarell", 0, 15)); // NOI18N
        EoutLabel.setForeground(java.awt.Color.blue);
        EoutLabel.setText("E out:");
        EoutLabel.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                EoutLabelMouseClicked(evt);
            }
        });

        EoutTextField.setEditable(false);
        EoutTextField.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseReleased(java.awt.event.MouseEvent evt) {
                EoutTextFieldMouseReleased(evt);
            }
        });
        EoutTextField.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                EoutTextFieldActionPerformed(evt);
            }
        });

        ThicknessLabel.setText("Thickness");
        ThicknessLabel.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                ThicknessLabelMouseClicked(evt);
            }
        });

        ThicknessTextField.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseReleased(java.awt.event.MouseEvent evt) {
                ThicknessTextFieldMouseReleased(evt);
            }
        });
        ThicknessTextField.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ThicknessTextFieldActionPerformed(evt);
            }
        });

        jLabel6.setText("um");

        jLabel7.setText("MeV");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(EinLabel)
                        .addGap(50, 50, 50)
                        .addComponent(EinTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 90, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(2, 2, 2)
                        .addComponent(jLabel2))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jLabel3)
                        .addGap(32, 32, 32)
                        .addComponent(modeComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(javax.swing.GroupLayout.Alignment.LEADING, layout.createSequentialGroup()
                                .addComponent(EoutLabel)
                                .addGap(36, 36, 36)
                                .addComponent(EoutTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 91, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addGroup(javax.swing.GroupLayout.Alignment.LEADING, layout.createSequentialGroup()
                                .addComponent(ThicknessLabel)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(ThicknessTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 90, javax.swing.GroupLayout.PREFERRED_SIZE)))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabel6)
                            .addComponent(jLabel7, javax.swing.GroupLayout.Alignment.TRAILING))))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(EinLabel)
                    .addComponent(EinTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel2))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(modeComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel3))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(ThicknessLabel)
                    .addComponent(ThicknessTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel6))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(EoutLabel)
                    .addComponent(EoutTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel7))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
    }// </editor-fold>//GEN-END:initComponents

    private void EinTextFieldActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_EinTextFieldActionPerformed
        exec.execute(this);
    }//GEN-LAST:event_EinTextFieldActionPerformed

    private void modeComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_modeComboBoxActionPerformed
        exec.execute(this);
        
        // also update the label for units:
        int modeIndex = modeComboBox.getSelectedIndex();
        String units = (String)modeComboBox.getItemAt(modeIndex);
        
        // remove MeV/ and parens:
        units = units.replace("MeV/","");
        units = units.replace("(","");
        units = units.replace(")","");
        units = units.trim();
        jLabel6.setText(units);
    }//GEN-LAST:event_modeComboBoxActionPerformed

    private void EoutTextFieldActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_EoutTextFieldActionPerformed
        exec.execute(this);
    }//GEN-LAST:event_EoutTextFieldActionPerformed

    private void ThicknessTextFieldActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ThicknessTextFieldActionPerformed
        exec.execute(this);
    }//GEN-LAST:event_ThicknessTextFieldActionPerformed

    private void EoutLabelMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_EoutLabelMouseClicked
        chooseMode(evt);
    }//GEN-LAST:event_EoutLabelMouseClicked

    private void jLabel3MouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_jLabel3MouseClicked
        // TODO add your handling code here:
    }//GEN-LAST:event_jLabel3MouseClicked

    private void ThicknessLabelMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_ThicknessLabelMouseClicked
        chooseMode(evt);
    }//GEN-LAST:event_ThicknessLabelMouseClicked

    private void EinLabelMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_EinLabelMouseClicked
        chooseMode(evt);
    }//GEN-LAST:event_EinLabelMouseClicked

    private void EinTextFieldMouseReleased(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_EinTextFieldMouseReleased
        // check for right mouse button and mode:
        if( evt.getButton() == java.awt.event.MouseEvent.BUTTON3 
                && calcMode == MODE_EIN )
        {
            // calculate coordinates:
            int x = evt.getX() + EinTextField.getX();
            int y = evt.getY() + EinTextField.getY();
            rightClickMenu.show(this, x, y);
        }
    }//GEN-LAST:event_EinTextFieldMouseReleased

    private void ThicknessTextFieldMouseReleased(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_ThicknessTextFieldMouseReleased
        // check for right mouse button and mode:
        if( evt.getButton() == java.awt.event.MouseEvent.BUTTON3 
                && calcMode == MODE_THICKNESS )
        {
            // calculate coordinates:
            int x = evt.getX() + ThicknessTextField.getX();
            int y = evt.getY() + ThicknessTextField.getY();
            rightClickMenu.show(this, x, y);
        }
    }//GEN-LAST:event_ThicknessTextFieldMouseReleased

    private void EoutTextFieldMouseReleased(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_EoutTextFieldMouseReleased
        // check for right mouse button and mode:
        if( evt.getButton() == java.awt.event.MouseEvent.BUTTON3 
                && calcMode == MODE_EOUT )
        {
            // calculate coordinates:
            int x = evt.getX() + EoutTextField.getX();
            int y = evt.getY() + EoutTextField.getY();
            rightClickMenu.show(this, x, y);
        }
    }//GEN-LAST:event_EoutTextFieldMouseReleased

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel EinLabel;
    private javax.swing.JTextField EinTextField;
    private javax.swing.JLabel EoutLabel;
    private javax.swing.JTextField EoutTextField;
    private javax.swing.JLabel ThicknessLabel;
    private javax.swing.JTextField ThicknessTextField;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel6;
    private javax.swing.JLabel jLabel7;
    private javax.swing.JComboBox modeComboBox;
    // End of variables declaration//GEN-END:variables

}
