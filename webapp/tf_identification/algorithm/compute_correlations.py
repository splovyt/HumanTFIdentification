import pandas as pd
from scipy.stats import spearmanr, pearsonr
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from time import sleep, time
import sys, os, pickle



_, GE_filepath, do_gorilla, amount_top_TFs = sys.argv
BASE_DIR = os.getcwd()

do_gorilla = bool(int(do_gorilla))
amount_top_TFs = int(amount_top_TFs)

# test sample:
# GE_filepath = 'aging_GE_data.csv'

def test_correlation(motif_df, expression_df, absolute_corr=False, absolute_expr=True):
    ''' Test the pairwise correlation between a df with motif z-score vectors after promoter scanning
    and a custom expression dataframe'''

    ## take the absolute values of the expression matrix (default=True)
    if absolute_expr:
        expression_df = abs(expression_df)

    ## Make rows equal
    index_intersection = set(motif_df.index) & set(expression_df.index)

    x = list()
    for idx in expression_df.index:
        if idx in index_intersection:
            x.append(True)
        else:
            x.append(False)
    expression_df = expression_df[x]
    expression_df['symbols'] = expression_df.index
    expression_df.drop_duplicates(subset='symbols', keep='first', inplace=True)
    expression_df.drop('symbols', inplace=True, axis=1)

    x = list()
    for idx in motif_df.index:
        if idx in index_intersection:
            x.append(True)
        else:
            x.append(False)
    motif_df = motif_df[x]

    ## test pairwise with spearman's rank correlation
    correlations = pd.DataFrame()
    for expr in expression_df.columns:
        c = list()
        for mtf in motif_df.columns:
            ## take the absolute values of the correlation vector (default=False)
            if absolute_corr:
                c.append(abs(pearsonr(motif_df[mtf], expression_df[expr])[0]))
            else:
                c.append(pearsonr(motif_df[mtf], expression_df[expr])[0])
        correlations[expr] = c
    correlations['Symbols'] = motif_df.columns
    correlations.set_index('Symbols', inplace=True)
    correlations.sort_values(correlations.columns[0], inplace=True, ascending=False)
    correlations.index = [x.replace('_zscore', '') for x in correlations.index]
    return correlations

def start_GOrilla(topTFs_set, allTFs_set):
    # Initialize the headless browser
    chrome_options = Options()
    chrome_options.add_argument("headless")
    try:
        os.system('chmod a+x ' + BASE_DIR + '/tf_identification/algorithm/chromedriver_linux')
        driver = webdriver.Chrome(BASE_DIR + '/tf_identification/algorithm/chromedriver_linux',
                              chrome_options=chrome_options)
    except: driver = webdriver.Chrome(BASE_DIR + '/tf_identification/algorithm/chromedriver_mac',
                              chrome_options=chrome_options)

    # Go to the GOrilla web application
    driver.get('http://cbl-gorilla.cs.technion.ac.il/')
    sleep(3)

    # Select the use of a background list
    driver.find_element_by_xpath('/html/body/form/blockquote[1]/table[2]/tbody/tr[2]/td/font[2]/input').click()

    # Select all GO categories
    driver.find_element_by_xpath('//*[@id="table4"]/tbody/tr[2]/td/input[4]').click()

    # Input top TFs
    topTFs = ''
    for tf in topTFs_set:
        topTFs += tf + '\n'
    element = driver.find_element_by_xpath('//*[@id="table3"]/tbody/tr[3]/td/textarea')
    element.clear()
    element.send_keys(topTFs)

    # Enter all TFs
    allTFs = ''
    for tf in allTFs_set:
        allTFs += tf + '\n'
    element = driver.find_element_by_xpath('//*[@id="table3"]/tbody/tr[7]/td/textarea')
    element.clear()
    element.send_keys(allTFs)

    sleep(15)

    # Submit query
    driver.find_element_by_xpath('/html/body/form/blockquote[2]/p/font/input').click()

    sleep(5)
    url = driver.current_url
    driver.quit()
    return url

def load_motif_vectors(filepath):
    ''' Load our own motif vectors. '''
    return pd.read_pickle(filepath)

def load_expression_FC(filepath):
    ''' Load the expression FC data. '''
    # try CSV
    # TODO: detect header and index
    try:
        data = pd.read_csv(filepath, header=0)
        data = data.set_index(data.columns[0])
        return data
    except: pass

    # try XLS
    try: return pd.read_excel(filepath)
    except: pass

    return pd.DataFrame()

def write_ERROR_html(cause):
    assert cause in ['input', 'calculation'], 'The ERROR cause should be in the list of known causes. Supply a correct cause.'
    # TODO: implement this
    # open the error_template.html
    with open(BASE_DIR + '/tf_identification/algorithm/' + 'error_template.html', 'r') as template_html:
        # add changes
        if cause == 'input':
            new_html = [line.replace('[TEXT GOES HERE]', 'The input file was not read correctly.') for line in template_html]
        elif cause == 'calculation':
            new_html = [line.replace('[TEXT GOES HERE]', 'Something went wrong.<br>Make sure your file matches the required format.') for line in
                        template_html]
    # save it in the same folder under the process_ID.html name
    with open(BASE_DIR + '/tf_identification/algorithm/' + str(process_ID) + '.html', 'w') as output_html:
        for line in new_html:
            output_html.write(line)
    # transfer it to the results folder for server access
    os.rename(BASE_DIR + '/tf_identification/algorithm/' + process_ID+'.html', BASE_DIR + '/tf_identification/templates/tf_identification/results/{}.html'.format(process_ID))
    return


def write_json(json_dict, nr):
    filename = BASE_DIR + '/tf_identification/static/tf_identification/results/barchart/{}.json'.format(process_ID + str(nr))
    static_name = '/static/tf_identification/results/barchart/{}.json'.format(process_ID + str(nr))
    with open(filename,
              'w') as outfile:
        outfile.write('[\n')
        for idx, key in enumerate(json_dict.keys()):
            if idx+1 == len(json_dict.keys()):
                outfile.write('\t{"tf": "'+key+'", "value": '+str(json_dict[key])+'}\n]')
            else:
                outfile.write('\t{"tf": "'+key+'", "value": '+str(json_dict[key])+'},\n')
    return static_name

def generate_accordion_div(correlation_matrix, top=amount_top_TFs, GOrilla = do_gorilla):

    # load in the symbol to name mappings
    with open(BASE_DIR + '/tf_identification/algorithm/symbol_to_name_dict.pickle', 'rb') as infile:
        symbol_to_name_dict = pickle.load(infile)

    complete_div = ''
    for condition_nr, condition in enumerate(correlation_matrix.columns):
        print('generating accordion for condition:', condition)
        # template div

        condition_template = \
            '''           <div class="card">
                    <div class="card-header" id="heading{}">
                      <h5 class="mb-0">
                        <button class="btn btn-link" data-toggle="collapse" data-target="#collapse{}" aria-expanded="true" aria-controls="collapse{}" >
                          [CONDITION NAME]
                        </button>
                        [GORILLA]
                      </h5>
                    </div>

                    <div id="collapse{}" class="collapse show" aria-labelledby="heading{}" data-parent="#accordion">
                      <div class="card-body">
                        [GRAPH]
                        [CONDITION CONTENT]
                      </div>
                    </div>
                  </div>
            '''.format(condition_nr, condition_nr, condition_nr, condition_nr, condition_nr)

        # sort TFs
        correlation_matrix.sort_values(condition, inplace=True, ascending=False)
        # extract top TFs
        topTFs = correlation_matrix.index[0:top]
        topCorr = correlation_matrix[condition][0:top]

        # if GOrilla is specified
        GOrilla_query_url = ''
        if GOrilla:
            topTFs_set = set([x.split('_')[0] for x in topTFs])
            allTFs_set = set([x.split('_')[0] for x in correlation_matrix.index])
            try: GOrilla_query_url = start_GOrilla(topTFs_set, allTFs_set)
            except: GOrilla_query_url = ''

        # convert to HTML
        topTFs_html = """<table class="table table-striped" style="text-align: center">
  <thead>
    <tr>
      <th scope="col">#</th>
      <th scope="col">Corr.</th>
      <th scope="col">HGNC</th>
      <th scope="col">Transcription factor name</th>
      <th scope="col">Motif</th>
      <th scope="col">WebLogo</th>
    </tr>
  </thead>
  <tbody style='font-size: 15px;'>
"""

        json_dict = {}
        for idx, tf in enumerate(topTFs):
            tf = tf.replace('_', ' ')
            try: tf_name = symbol_to_name_dict[tf.split(' ')[0]]
            except: tf_name = 'NA'
            # add motif logo
            TF_img = '<img src="/static/tf_identification/img/motiflogos/{}.png" width="100px">'.format(tf.replace(' ', '_'))
            # add features in table format
            topTFs_html += """
                 <tr>
      <th scope="row">{}</th>
      <td>{}</td>
      <td>{}</td>
      <td>{}</td>
      <td>{}</td>
      <td>{}</td>
    </tr>""".format(idx+1, round(topCorr[idx],3), tf.split(' ')[0], tf_name, tf.split(' ')[1], TF_img)
            json_dict[tf] = round(topCorr[idx],3)
        topTFs_html += """
          </tbody>
</table>"""

        # produce chart
        # write json for chart
        json_static_filename = write_json(json_dict, condition_nr)
        # write html for graph
        graph_html = ''
        with open(BASE_DIR + '/tf_identification/static/tf_identification/results/barchart/chart_template.html') as infile:
            for line in infile:
                graph_html += line.replace('[JSON FILE]', json_static_filename).replace('[NR]', str(condition_nr))


        # insert in template
        # TODO: add the condition_nr to auto-fold everything but the first one
        if GOrilla_query_url != '':
            GOrilla_query_url = "<a class='btn btn-primary float-right' href='{}' target='_blank'>Open GOrilla results</a>".format(GOrilla_query_url)

        complete_div += condition_template.replace('[CONDITION NAME]', condition) \
            .replace('[CONDITION CONTENT]', topTFs_html) \
            .replace('[GRAPH]', graph_html).replace('class="collapse show"', 'class="collapse"')\
            .replace('[GORILLA]', GOrilla_query_url)

    return complete_div

def write_SUCCESS_html(correlations):
    # TODO: implement this
    correlations.to_csv(
        BASE_DIR + '/tf_identification/static/tf_identification/results/{}.csv'.format(process_ID),
        header=True,
        index=True)
    # open the results_template.html
    with open(BASE_DIR + '/tf_identification/algorithm/' + 'results_template.html', 'r') as template_html:
        # add changes
        condition = correlations.columns[0]
        top_TF = condition + '<br></br>'
        for TF in correlations.index[0:5]:
            top_TF += TF + '<br>'

        new_html = []
        for line in template_html:
            # download link
            line = line.replace('[DOWNLOAD LINK]', '/static/tf_identification/results/{}.csv'.format(process_ID))
            # replace TOP TFS in here
            if '[ACCORDION]' in line:
                line = line.replace('[ACCORDION]', generate_accordion_div(correlations))
            # append
            new_html.append(line)

    # save it in the same folder under the process_ID.html name
    with open(BASE_DIR + '/tf_identification/algorithm/' + str(process_ID) + '.html', 'w') as output_html:
        for line in new_html:
            output_html.write(line)
    # transfer it to the results folder for server access
    os.rename(BASE_DIR + '/tf_identification/algorithm/' + process_ID+'.html', BASE_DIR + '/tf_identification/templates/tf_identification/results/{}.html'.format(process_ID))
    return

def clean(directory):
    for filename in os.listdir(directory):
        # remove the file extension
        file = filename[0:filename.find('.')]

        # if the length of the file is equal to the length of the process ID
        if len(file) >= 29:
            try:
                file = file[0:30]
                # check creation time
                creation_time = int(file[-10:])
                # delete the file if needed
                if (time() - creation_time) >= storage_time:
                    os.remove(directory + filename)
            except: pass

# Load data
#motif_vectors = load_motif_vectors(BASE_DIR + '/tf_identification/algorithm/Li_lab_implementation_vectors.pickle')
motif_vectors = load_motif_vectors(BASE_DIR + '/tf_identification/algorithm/FIMO_log2.pickle')

process_ID = GE_filepath[-30:]
expression_vectors = load_expression_FC(GE_filepath)
if expression_vectors.shape[0] == 0:
    write_ERROR_html(cause='input')

# Compute correlation
try:
    correlation_matrix = test_correlation(motif_vectors, expression_vectors,
                                      absolute_corr=False, absolute_expr=True)
    write_SUCCESS_html(correlation_matrix)
except:
    write_ERROR_html(cause='calculation')

# Clean up X days old files
storage_time = 60*60*24*31 # 1 month in seconds

# directories to clean up:
# /tf_identification/static/tf_identification/results/
clean(BASE_DIR + '/tf_identification/static/tf_identification/results/')
# /tf_identification/static/tf_identification/results/barchart/
clean(BASE_DIR + '/tf_identification/static/tf_identification/results/barchart/')
# /tf_identification/templates/tf_identification/results/
clean(BASE_DIR + '/tf_identification/templates/tf_identification/results/')
# /useruploads/
clean(BASE_DIR + '/useruploads/')