import QtQuick 2.8
import JASP.Controls 1.0
import JASP.Theme 1.0
import JASP.Widgets 1.0

Form
{
  	VariablesForm
  	{
  		height: 300
  		AvailableVariablesList { name: "allVariablesList" }		
  		AssignedVariablesList { name: "variables"; title: qsTr("Variables")}
  	}
  	
  	Section
    {
      title: qsTr("Single-Test Reliability")
    	Group
    	{
    		title: qsTr("Scale Statistics")
    		CheckBox 
    		{    
    		  name:   "mcDonaldScale"	
    		  label:  qsTr("McDonald's ω")         
    		  id:     mcdonald
    
    		  CheckBox 
    		  {    
    		  name:   "dispPPC"	
    		  label:  qsTr("Posterior predictive check"); 
    		  enabled: mcdonald.checked        
    	  	}
    
    		}
        CheckBox 
        {     
          name: "alphaScale";				
          label: qsTr("Cronbach's α");         
          id: cronbach   		      
          }
        CheckBox 
        {     
          name: "guttman2Scale";			
          label: qsTr("Guttman's λ2");        
          id: guttman2      
          }
 //       CheckBox 
   //     {     
   //       name: "guttman6Scale";			
   //       label: qsTr("Guttman's λ6");        
   //       id: guttman6      
  //        }
    		CheckBox 
    		{     
    		  name: "glbScale";				  
    		  label: qsTr("Greatest lower bound"); 
    		  id: glb      	              
    		  }
	      CheckBox { name: "averageInterItemCor";				label: qsTr("Average interitem correlation")	}

    		CIField 
    		{      
    		name: "CredibleIntervalValue";   
    		label: qsTr("Credible interval");
    		defaultValue: 95
    		}
        CheckBox { name: "meanScale";						label: qsTr("Mean")							}
        CheckBox { name: "sdScale";							label: qsTr("Standard deviation")			}
        
    	}
           
    	Group
    	{
    		title: qsTr("Individual Item Statistics")
    		CheckBox 
    		{ 
    		  name: "mcDonaldItem";				
    		  label: qsTr("McDonald's ω  (if item dropped)");	        
    		  enabled: mcdonald.checked 
    		  }
    		CheckBox 
    		{ 
    		  name: "alphaItem";					
    		  label: qsTr("Cronbach's α (if item dropped)");	        
    		  enabled: cronbach.checked 
    		  }
    		CheckBox 
    		{ 
    		  name: "guttman2Item";				
    		  label: qsTr("Guttman's λ2 (if item dropped)");	        
    		  enabled: guttman2.checked  
    		  }
//    		CheckBox 
//    		{ 
//    		  name: "guttman6Item";				
//    		  label: qsTr("Guttman's λ6 (if item dropped)");	        
//    		  enabled: guttman6.checked  
//    		  }
        CheckBox 
        { 
          name: "glbItem";     				
          label: qsTr("Greatest lower bound (if item dropped)");	
          enabled: glb.checked     
          }
        CheckBox { name: "itemRestCor";						label: qsTr("Item-rest correlation")				}

        CheckBox { name: "meanItem";						label: qsTr("Mean")								}
        CheckBox { name: "sdItem";							label: qsTr("Standard deviation")				}
    	}
    
      Group
      {
          CheckBox 
          {
            name: "plotPosterior";           
            label: qsTr("Plot Posteriors");
            CheckBox 
            {   
              name: "fixXRange";               
              label: qsTr("Fix range to 0-1")
              }
            CheckBox { 
              name: "dispPrior";               
              label: qsTr("Display Priors")
              }
        
            }
        }
        
        Group
        {
            CheckBox 
            { 
              id:                 probTable
              name:               "probTable"; 
              label:              qsTr("Probability for: ")
              childrenOnSameRow:  true
              
              DoubleField
              {
                name: "probTableValueLow"
                afterLabel: qsTr("< Reliability <")
                defaultValue: 0.70
                min: 0
                max: 1
                decimals: 2
                fieldWidth: 40
              }
              DoubleField
              {
                name: "probTableValueHigh"
               // label: qsTr("< Reliability <")
                defaultValue: 1
                min: 0
                max: 1
                decimals: 2
                fieldWidth: 40
                  
              }
              
            }
            Item
            {
              width:  shadePlots.width + Theme.subOptionOffset
              height: shadePlots.height
              
              CheckBox 
              { 
                indent:   true
                id:       shadePlots
                name:     "shadePlots";              
                label:    qsTr("Shade posterior region in plots"); 
                enabled:  probTable.checked    
                x:        Theme.subOptionOffset
              }
            }
          }
      }
      
  
	Section
	{
		title: qsTr("Reverse-Scaled Items")
		
		VariablesForm
		{
			height: 150
			AvailableVariablesList { name: "normalScaledItems";	 title: qsTr("Normal-Scaled Items"); source: "variables" }
			AssignedVariablesList {  name: "reverseScaledItems"; title: qsTr("Reverse-Scaled Items") }
		}
	}
	
	Section
	{
		title: qsTr("Advanced Options")
        IntegerField
        {
            name: "noSamples"
            label: qsTr("No. of posterior samples")
            defaultValue: 500
            fieldWidth: 50
            min: 100
            max: 1e7
        }
//        RadioButtonGroup {
//            title: qsTr("Missing Values")
//            name: "missingValues"
//            RadioButton { value: "excludeCasesListwise"; label: qsTr("Exclude cases listwise"); checked: true	}
//            RadioButton { value: "excludeCasesPairwise"; label: qsTr("Exclude cases pairwise")					}
//        }

    }
}
