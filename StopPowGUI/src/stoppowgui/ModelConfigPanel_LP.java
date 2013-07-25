package stoppowgui;

import stoppowgui.util.TableFocusAction;
import java.awt.event.KeyEvent;
import javax.swing.table.DefaultTableModel;
import cStopPow.StopPow_LP;
import cStopPow.FloatVector;
import SciTK.DialogError;
import java.awt.event.ActionEvent;
import javax.swing.AbstractAction;
import javax.swing.KeyStroke;


/**
 * Implement a dialog to prompt the user for values
 * to construct a new Li-Petrasso stopping power model.
 * The model is automatically added to the model
 * manager at the end.
 * 
 * @brief Initialize a L-P model
 * @class ModelConfigPanel_LP
 * @author Alex Zylstra
 * @date 2013/06/05
 */
public class ModelConfigPanel_LP extends ModelConfigPanel {
    /**
     * Creates new form ModelConfigPanel_LP
     * @param parent the JFrame parent of this dialog
     * @param modal the modal mode for this dialog
     * @param m the ModelManager in use by this application (i.e. where the new Li-Petrasso model should be added)
     */
    public ModelConfigPanel_LP(ModelManager m) {
        super(m);
        
        initComponents();
        TableFocusAction myTableFocusAction = new TableFocusAction(table_plasma, KeyStroke.getKeyStroke("TAB"), this.text_name);
        
        // add some key bindings:
        button_cancel.getInputMap().put(KeyStroke.getKeyStroke("released ENTER"), "Enter");
        button_cancel.getActionMap().put("Enter", new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e)
            {
                button_cancelActionPerformed(e);
            }
        });
        button_ok.getInputMap().put(KeyStroke.getKeyStroke("released ENTER"), "Enter");
        button_ok.getActionMap().put("Enter", new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e)
            {
                button_okActionPerformed(e);
            }
        });
        
        setVisible(true);
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jScrollPane1 = new javax.swing.JScrollPane();
        table_plasma = new javax.swing.JTable();
        jPanel1 = new javax.swing.JPanel();
        jLabel1 = new javax.swing.JLabel();
        jLabel2 = new javax.swing.JLabel();
        jLabel3 = new javax.swing.JLabel();
        jLabel4 = new javax.swing.JLabel();
        text_name = new javax.swing.JTextField();
        text_A = new javax.swing.JTextField();
        text_Z = new javax.swing.JTextField();
        checkbox_collective = new javax.swing.JCheckBox();
        button_ok = new javax.swing.JButton();
        button_cancel = new javax.swing.JButton();

        setPreferredSize(new java.awt.Dimension(370, 350));
        setSize(new java.awt.Dimension(370, 350));
        setLayout(new org.netbeans.lib.awtextra.AbsoluteLayout());

        jScrollPane1.setPreferredSize(new java.awt.Dimension(350, 200));

        table_plasma.setModel(new javax.swing.table.DefaultTableModel(
            new Object [][] {
                {null, null, null, null},
                {null, null, null, null}
            },
            new String [] {
                "A", "Z", "n [1/cc]", "T [keV]"
            }
        ) {
            Class[] types = new Class [] {
                java.lang.Float.class, java.lang.Float.class, java.lang.Float.class, java.lang.Float.class
            };

            public Class getColumnClass(int columnIndex) {
                return types [columnIndex];
            }
        });
        table_plasma.setColumnSelectionAllowed(true);
        table_plasma.setGridColor(new java.awt.Color(0, 0, 0));
        table_plasma.setNextFocusableComponent(button_ok);
        table_plasma.getTableHeader().setReorderingAllowed(false);
        table_plasma.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyPressed(java.awt.event.KeyEvent evt) {
                table_plasmaKeyPressed(evt);
            }
        });
        jScrollPane1.setViewportView(table_plasma);
        table_plasma.getColumnModel().getSelectionModel().setSelectionMode(javax.swing.ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        table_plasma.getColumnModel().getColumn(0).setResizable(false);
        table_plasma.getColumnModel().getColumn(1).setResizable(false);
        table_plasma.getColumnModel().getColumn(2).setResizable(false);
        table_plasma.getColumnModel().getColumn(3).setResizable(false);

        add(jScrollPane1, new org.netbeans.lib.awtextra.AbsoluteConstraints(10, 5, -1, -1));

        jPanel1.setPreferredSize(new java.awt.Dimension(350, 130));

        jLabel1.setText("A:");

        jLabel2.setText("Z:");

        jLabel3.setText("Test particle");

        jLabel4.setText("Name:");

        text_name.setText("Li-Petrasso");
        text_name.setNextFocusableComponent(text_A);

        checkbox_collective.setText("Collective effects");

        button_ok.setText("OK");
        button_ok.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                button_okActionPerformed(evt);
            }
        });

        button_cancel.setText("Cancel");
        button_cancel.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                button_cancelActionPerformed(evt);
            }
        });

        org.jdesktop.layout.GroupLayout jPanel1Layout = new org.jdesktop.layout.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(jPanel1Layout.createSequentialGroup()
                .add(jPanel1Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(jPanel1Layout.createSequentialGroup()
                        .addContainerGap()
                        .add(jLabel3)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(jLabel1)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.UNRELATED)
                        .add(jPanel1Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                            .add(checkbox_collective)
                            .add(jPanel1Layout.createSequentialGroup()
                                .add(text_A, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 84, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                .add(jLabel2)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                .add(text_Z, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 84, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))))
                    .add(jPanel1Layout.createSequentialGroup()
                        .add(94, 94, 94)
                        .add(button_ok)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(button_cancel))
                    .add(jPanel1Layout.createSequentialGroup()
                        .addContainerGap()
                        .add(jLabel4)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.UNRELATED)
                        .add(text_name, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 141, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(43, Short.MAX_VALUE))
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(jPanel1Layout.createSequentialGroup()
                .add(0, 0, 0)
                .add(jPanel1Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel4)
                    .add(text_name, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(jPanel1Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(jLabel1)
                    .add(text_A, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 28, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(jLabel2)
                    .add(text_Z, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 28, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(jLabel3))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(checkbox_collective)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED, 10, Short.MAX_VALUE)
                .add(jPanel1Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(button_ok)
                    .add(button_cancel)))
        );

        add(jPanel1, new org.netbeans.lib.awtextra.AbsoluteConstraints(10, 211, -1, -1));
    }// </editor-fold>//GEN-END:initComponents

    /**
     * Handle dynamic allocation of additional rows in the table
     * via key events. Also, handles tab navigation from the table
     * to the next text component.
     * @param evt The KeyEvent that triggered this call
     */
    private void table_plasmaKeyPressed(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_table_plasmaKeyPressed
        // info on where we are in the table:
        int num_rows = table_plasma.getRowCount();
        int num_cols = table_plasma.getColumnCount();
        int row = table_plasma.getSelectedRow();
        int col = table_plasma.getSelectedColumn();
        
        // whether we should add a new row:
        boolean add_row = false;
        
        // first, check if user typed enter:
        if ( evt.getKeyCode() == KeyEvent.VK_ENTER )
        {
            // only trigger based on enter if we are at the lower
            // right hand corner of table
            // (do this because of how enter traverses JTables)
            if( row+1 == num_rows && col+1 == num_cols)
                add_row = true;
        } 
        // next, check for down arrow events:
        if (evt.getKeyCode() == KeyEvent.VK_DOWN )
        {
            // only trigger if we are at the bottom row:
            if( row+1 == num_rows )
                add_row = true;
        }
        // if necessary, add the row:
        if( add_row )
        {
            // add a new row to the bottom of the table:
            Float[] new_row = new Float[] {null,null,null,null};
            ( (DefaultTableModel) table_plasma.getModel() ).addRow(new_row);
        }
        
        // also check for tab
        if( evt.getKeyCode() == KeyEvent.VK_TAB && !table_plasma.isEditing() )
        {
            // if we are the lower right hand corner of the table,
            // help out the user by moving focus
            if( row+1 == num_rows && col+1 == num_cols )
            {
                // move cursor to the first text box:
                //this.text_name.requestFocusInWindow();
                // deselect cell of table:
                //table_plasma.getSelectionModel().clearSelection();
            }
            // otherwise, advance to the next cell:
            else
            {
                //Component next = table_plasma.getNextFocusableComponent();
                //next.requestFocus();
            }
        }
    }//GEN-LAST:event_table_plasmaKeyPressed

    /**
     * Action taken when the "OK" button is clicked. What this does
     * is validate the entered parameters, and generate an error
     * if they are not OK.
     * If this dialog was created from a ModelManagerGUI, then the model
     * is automatically created and added to the parent.
     * @param evt The ActionEvent triggering this call
     */
    private void button_okActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_button_okActionPerformed
        // first, validate results:
        boolean data_ok = validate_data();
        if( !data_ok )
        {
            DialogError p = new DialogError(this,"Values entered are not OK.");
            return;
        }
        
        // if we can, update the ModelManager automatically:
        if( models != null )
        {
            // get the model:
            StopPow_LP s = get_model();
            
            // get the name:
            String name = text_name.getText();
            
            // if we are not in edit mode, and the name is unique:
            if( !models.containsKey(name) )
            {
                // Add new model to the ModelManagerGUI:
                models.add_model(name, s, this);
                previousKey = name;

            }
            else if( editMode && name.equals(previousKey) )
            {
                models.change_model(name, s, this);
            }
            else
            { // name is not unique
                DialogError p = new DialogError(this,"Name is not unique");
                return;
            }
            
            // set this panel to edit mode:
            editMode = true;
            
        }
        dialog.dispose();
    }//GEN-LAST:event_button_okActionPerformed

    /**
     * Action taken when the "cancel" button is clicked, i.e. dispose
     * of this window.
     * @param evt The ActionEvent triggering this call
     */
    private void button_cancelActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_button_cancelActionPerformed
        // get rid of the window:
        dialog.dispose();
    }//GEN-LAST:event_button_cancelActionPerformed

    // ---------------------------------------
    //  Methods for dealing with entered data
    // ---------------------------------------
    /**
     * Validate the data contained within the dialog.
     * @return true if all entered parameters are OK, false otherwise.
     */
    public boolean validate_data()
    {
        // get access to table data:
        DefaultTableModel table_data = (DefaultTableModel) table_plasma.getModel();
        
        try
        {
            // iterate over rows of the table:
            for(int i=0; i < table_data.getRowCount(); i++)
            {
                // check A value:
                Float A = (Float)table_data.getValueAt(i, 0);
                if( A == null || A <= 0 )
                    return false;

                // Check Z value:
                Float Z = (Float)table_data.getValueAt(i, 1);
                if( Z == null || Z <= 0 )
                    return false;

                // Check n value:
                Float n = (Float)table_data.getValueAt(i, 2);
                if( n == null || n <= 0 )
                    return false;

                // Check Z value:
                Float T = (Float)table_data.getValueAt(i, 3);
                if( T == null || T <= 0 )
                    return false;
            }

            // check the name:
            String name = text_name.getText();
            if( name.equals("") )
                return false;

            // check the test particle Z
            Float Zt = Float.parseFloat( text_Z.getText() );
            if( Zt == 0 )
                return false;

            // check the test particle A
            Float At = Float.parseFloat( text_A.getText() );
            if( At <= 0 )
                return false;
        }
        catch(Exception e)
        {
            return false;
        }
        
        return true;
    }
    
    /**
     * Construct a new Li-Petrasso model reflecting the values entered.
     * @return new stopping power model
     */
    public StopPow_LP get_model()
    {
        // This method will construct a new Li-Petrasso stopping power model
        
        // get access to table data:
        DefaultTableModel table_data = (DefaultTableModel) table_plasma.getModel();
        
        // Get number of field particles:
        int num = table_data.getRowCount();
        
        // Construct FloatVectors for making the StopPow_LP:
        FloatVector Af = new FloatVector(num);
        FloatVector Zf = new FloatVector(num);
        FloatVector nf = new FloatVector(num);
        FloatVector Tf = new FloatVector(num);
        
        // loop over all rows of the table to populate the above:
        for(int i=0; i < num; i++)
        {
            Af.set( i , (Float)table_data.getValueAt(i, 0) );
            Zf.set( i , (Float)table_data.getValueAt(i, 1) );
            nf.set( i , (Float)table_data.getValueAt(i, 2) );
            Tf.set( i , (Float)table_data.getValueAt(i, 3) );
        }
        
        // get the test particle info:
        Float Zt = Float.parseFloat( text_Z.getText() );
        Float At = Float.parseFloat( text_A.getText() );
        
        // create the StopPow_LP object:
        StopPow_LP model;
        try
        {
            model = new StopPow_LP(At,Zt,Af,Zf,Tf,nf);
        }
        catch(java.lang.IllegalArgumentException e)
        {
            DialogError d = new DialogError(this,"Could not create model: " + e.getMessage());
            return null;
        }
        
        // Get the collective effects flag, and set it:
        boolean collective = checkbox_collective.isSelected();
        model.set_collective(collective);
        
        return model;
    }
    
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton button_cancel;
    private javax.swing.JButton button_ok;
    private javax.swing.JCheckBox checkbox_collective;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JTable table_plasma;
    private javax.swing.JTextField text_A;
    private javax.swing.JTextField text_Z;
    private javax.swing.JTextField text_name;
    // End of variables declaration//GEN-END:variables
    private ModelManagerGUI model_manager;

    @Override
    public String get_type() {
        return "Li-Petrasso";
    }
}
