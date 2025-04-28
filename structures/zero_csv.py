import pandas as pd

def create_zeroed_csv(output_path,rows):
    columns = ['Deflection', 'Twist']
    row_count = rows

    zero_df = pd.DataFrame(0, index=range(row_count), columns=columns)

    # Write to CSV
    zero_df.to_csv(output_path, index=False)

rows = 11
create_zeroed_csv('def_twist_pos.csv',rows)
create_zeroed_csv('def_twist_neg.csv',rows)