import dash
from dash import dcc, html, Input, Output, State, dash_table
import pandas as pd
import base64
import io
import os

# Initialize Dash app
app = dash.Dash(__name__)
server = app.server  # Needed for deployment

# Layout of the web app
app.layout = html.Div([
    html.H1("FastQ File Processor"),
    dcc.Upload(
        id='upload-data',
        children=html.Button('Upload File'),
        multiple=False
    ),
    html.Div(id='output-data-upload'),
])


# Callback to process uploaded file
@app.callback(
    Output('output-data-upload', 'children'),
    Input('upload-data', 'contents'),
    State('upload-data', 'filename')
)
def update_output(contents, filename):
    if contents is None:
        return "No file uploaded."

    # Decode file
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)

    # Save file temporarily
    temp_filename = f"temp_{filename}"
    with open(temp_filename, "wb") as f:
        f.write(decoded)

    # Process file (placeholder for real processing)
    result = pd.DataFrame({"Filename": [filename], "Status": ["Processed"]})

    # Remove temp file after processing
    os.remove(temp_filename)

    return dash_table.DataTable(
        data=result.to_dict('records'),
        columns=[{"name": i, "id": i} for i in result.columns]
    )


# Run app
if __name__ == '__main__':
    app.run_server(debug=True)
