

public class main {
	public static void main(String argv[]) {
		System.loadLibrary("cStopPow"); // load the JNI

		// various models:
		StopPow[] models = new StopPow[3];

		// set up SRIM calculator:
		// protons in normal solid aluminum
		String fname = new String("data/Hydrogen in Aluminum.txt");
		StopPow s = new StopPow_SRIM(fname);
		models[0] = s;

		// set up Li-Petrasso
		// proton in hydrogen plasma at 1e23 density at 1keV temperature
		FloatVector mf = new FloatVector(2);
		mf.set( 0 , 1.0f );
		mf.set( 1 , 1/1800.f);
		FloatVector Zf = new FloatVector(2);
		Zf.set( 0 , 1 );
		Zf.set( 1 , -1);
		FloatVector Tf = new FloatVector(2);
		Tf.set( 0 , 1 );
		Tf.set( 1 , 1 );
		FloatVector nf = new FloatVector(2);
		nf.set( 0 , 1e23f );
		nf.set( 1 , 1e23f );
		StopPow s2 = new StopPow_LP(1,1,mf,Zf,Tf,nf);
		models[1] = s2;

		// set up Bethe-Bloch
		// proton in solid C (diamond)
		FloatVector mf2 = new FloatVector(1);
		mf2.set( 0 , 12.0f );
		FloatVector Zf2 = new FloatVector(1);
		Zf2.set( 0 , 6.0f );
		FloatVector nf2 = new FloatVector(1);
		nf2.set( 0 , 1.76e23f );
		StopPow s3 = new StopPow_BetheBloch(1,1,mf2,Zf2,nf2);
		models[2] = s3;

		// perform some tests
		for(int i=0; i < models.length; i++)
		{
			models[i].set_mode( models[i].getMODE_LENGTH() );
			System.out.println( "dEdx(10 MeV) = " + models[i].dEdx(10) + " MeV/um");
			models[i].set_mode( models[i].getMODE_RHOR() );
			System.out.println( "dEdx(10 MeV) = " + models[i].dEdx(10) + " MeV/(mg/cm2)");

			models[i].set_mode( models[i].getMODE_LENGTH() );
			System.out.println( "Eout(10 MeV, 100um) = " + models[i].Eout(10,100) );
			System.out.println( "Ein(10 MeV, 100um) = " + models[i].Ein(10,100) );
			System.out.println( "Thickness(10 MeV, 9 MeV) = " + models[i].Thickness(10,9) );

			System.out.println("--------------------");
		}
	}
}