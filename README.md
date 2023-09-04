# matplotlib-challenge
Matplotlib for pharmaceuticals

Analysis
One observation that can be made from the data is capomulin and ramicane drug performance are very similar. The summary statistics table depicts values for mean, median, variance, standard deviation, and standard error tumor volumes that don't deviate significantly between the two drug regimens. It can be inferred that capomulin and ramicane are the most recommended drug regimens to be used for treatment of mice with squamous cell carcinoma. Another observation is that the higher the mouse weight, the higher the tumor volume. This supports the idea that regardless of the drug's effectiveness, tumor development is still expected to increase for bigger-sized mice. A third observation is that ketapril has the highest mean tumor volume, which indicates that this drug is the least effective in treating mice with SCC.


#Specify imports for data visualization, pandas, and statistics
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st

#Sets the file path from the current, which is specified as ".", to the specific file
mouse_metadata_path = "./data/Mouse_metadata.csv"
study_results_path = "./data/Study_results.csv"

#Reads the mouse metadata and study results csv files through pandas
mouse_metadata_read = pd.read_csv(mouse_metadata_path)
study_results_read = pd.read_csv(study_results_path)

#Merge the datasets by using Mouse ID as the common denominator and includes the rows from the second dataset with the same mouse id and discards the other rows.
study_data_complete_df = pd.merge(study_results_read, mouse_metadata_read, how="left", on="Mouse ID")

#Show the merged dataset
study_data_complete_df


# In[ ]:


#Determines the number of unique mouse IDs in the dataset
len(study_data_complete_df["Mouse ID"].unique())


# In[ ]:


#Grab the duplicates of the dataset where mice id's have same timepoints. 
mouse_ID_duplicates = study_data_complete_df[study_data_complete_df.duplicated(subset=["Mouse ID", "Timepoint"])]["Mouse ID"].unique()
mouse_ID_duplicates


# In[ ]:


#Show data for the g989 duplicate mouse ID
dataset_mouse_duplicate = study_data_complete_df[study_data_complete_df["Mouse ID"] == "g989"]
dataset_mouse_duplicate


# In[ ]:


#Make a new dataframe where the g989 duplicate mouse ID is omitted
clean_data_df = study_data_complete_df[study_data_complete_df["Mouse ID"].isin(mouse_ID_duplicates) == False]
clean_data_df


# In[ ]:


#Determines the number of unique mouse IDs in the new clean dataset
len(clean_data_df["Mouse ID"].unique())


# ## Summary Statistics

# In[ ]:


#Provide summary statistics on the clean data dataframe
mean = clean_data_df.groupby("Drug Regimen").mean()["Tumor Volume (mm3)"]
median = clean_data_df.groupby("Drug Regimen").median()["Tumor Volume (mm3)"]
variance = clean_data_df.groupby("Drug Regimen").var()["Tumor Volume (mm3)"]
standard_deviation = clean_data_df.groupby("Drug Regimen").std()["Tumor Volume (mm3)"]
SEM = clean_data_df.groupby("Drug Regimen").sem()["Tumor Volume (mm3)"]

#Take the previous results, which are data series, and put them in a pandas dataframe.
summary_statistics_df = pd.DataFrame({
    "Mean Tumor Volume" : mean,
    "Median Tumor Volume" : median,
    "Tumor Volume Variance" : variance,
    "Tumor Volume Standard Deviation" : standard_deviation,
    "Tumor Volume Standard Error" : SEM,
})
#Display the results
summary_statistics_df


# In[ ]:


# A more advanced method to generate a summary statistics table of mean, median, variance, standard deviation,
# and SEM of the tumor volume for each regimen (only one method is required in the solution)

# Using the aggregation method, produce the same summary statistics in a single line


# ## Bar and Pie Charts

# In[ ]:


#Create bar graph for each drug regimen and total number of tested mice
counts_drug_reg = clean_data_df["Drug Regimen"].value_counts()
counts_drug_reg.plot(kind="bar")
plt.xlabel("Drug Regimen")
plt.ylabel("Number of Tested Mice")
plt.show()


# In[ ]:


#Create bar graph using pyplot
counts_drug_reg = clean_data_df["Drug Regimen"].value_counts()
plt.bar(counts_drug_reg.index.values, counts_drug_reg.values) 
plt.xlabel("Drug Regimen")
plt.xticks(rotation=90)
plt.ylabel("Number of Tested Mice")
plt.show()


# In[ ]:


#Pie chart using Pandas
count_drug_reg_pie = clean_data_df.Sex.value_counts()
count_drug_reg_pie.plot(kind="pie", autopct = "%1.1f%%")


# In[ ]:


#Same pie chart using pyplot
count_drug_reg_pie = clean_data_df.Sex.value_counts()
plt.pie(count_drug_reg_pie.values, labels=count_drug_reg_pie.index.values, autopct='%1.1f%%')
plt.ylabel("Sex")


# ## Quartiles, Outliers and Boxplots

# In[ ]:


#Calculate the greatest timepoint
greatest_tp = clean_data_df.groupby(["Mouse ID"])["Timepoint"].max()
greatest_tp = greatest_tp.reset_index()
# Merge this group df with the original DataFrame to get the tumor volume at the last timepoint
data_merged = greatest_tp.merge(clean_data_df, on=["Mouse ID", "Timepoint"], how="left")


# In[ ]:


#Create a list treatments with four specific drug regimens
treatments = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]

#Create empty list for tumor volume later on
tumor_vol = []

#For each drug in treatments list
for drug in treatments:
    
    #Determine the specific volume for each drug
    final_tvol = data_merged.loc[data_merged["Drug Regimen"] == drug, "Tumor Volume (mm3)"]
    
    #Add it to the empty list 
    tumor_vol.append(final_tvol)
    
    #Calculate outliers
    quartiles = final_tvol.quantile([0.25, .5, .75])
    lowerq = quartiles[0.25]
    upperq = quartiles[0.75]
    iqr = upperq - lowerq
    lowerbound = lowerq - (1.5* iqr)
    upperbound = upperq + (1.5 * iqr)
    outliers = final_tvol.loc[(final_tvol > lowerbound) | (final_tvol > upperbound)]


# In[ ]:


#Create a box plot showing the distribution of tumor volume for each treatment group
orange_out = dict(markerfacecolor='red', markersize=15)
plt.boxplot(tumor_vol, labels = treatments, flierprops=orange_out)
plt.ylabel("Final Tumor Volume (mm3)")
plt.show()


# ## Line and Scatter Plots

# In[ ]:


#Make a line for r554 mouse with capomulin drug regimen as an example
capomulin_data = clean_data_df[clean_data_df["Drug Regimen"] == "Capomulin"]
mouse_data = capomulin_data[capomulin_data["Mouse ID"] == "r554"]
plt.plot(mouse_data["Timepoint"], mouse_data["Tumor Volume (mm3)"])
plt.xlabel("Timepoint (days)")
plt.ylabel("Tumor Volume (mm3)")
plt.title("Capomulin treamtment of mouse r554 over time")


# In[ ]:


#Create a scatter plot of mouse weight and the average observed tumor volume for all those treated with Capomulin
capomulin_data = clean_data_df[clean_data_df["Drug Regimen"] == "Capomulin"]
capomulin_avg = capomulin_data.groupby(["Mouse ID"]).mean()
plt.scatter(capomulin_avg["Weight (g)"], capomulin_avg["Tumor Volume (mm3)"])
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.show()


# ## Correlation and Regression

# In[ ]:


#Find the correlation coefficient and a linear regression model for mouse weight and average observed tumor volume for all those trated with Capomulin
correlation = st.pearsonr(capomulin_avg["Weight (g)"], capomulin_avg["Tumor Volume (mm3)"])
print(f"The correlation between is {round(correlation[0],2)}")

#Create the x and y values
model = st.linregress(capomulin_avg["Weight (g)"], capomulin_avg["Tumor Volume (mm3)"])

#Set line slope to 0th index since the model displays slope as the first item in the list
slope = model[0]

#Set line slope to 1st index since the model displays slope as the first item in the list
b = model[1]
y_values = capomulin_avg["Weight (g)"] * slope + b
plt.scatter(capomulin_avg["Weight (g)"], capomulin_avg["Tumor Volume (mm3)"])
plt.plot(capomulin_avg["Weight (g)"], y_values, color="red")
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.show()

