from django.shortcuts import render, redirect
from django.views import View
from tf_identification.forms import UploadFileForm
from webapp.settings import BASE_DIR
from django.utils.crypto import get_random_string
import time, os, random

def handle_uploaded_file(f, filename, gorilla, topTF):
    ''' Handle and process the uploaded file accordingly. '''
    # Save the file from memory onto the server.
    target_file = BASE_DIR + '/useruploads/' + filename
    print('Uploaded file saved at', target_file)
    with open(target_file, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)
    # Launch the Python script on the file
    algorithm_script = BASE_DIR + '/tf_identification/algorithm/compute_correlations.py'
    os.system('python3 {} {} {} {} &'.format(algorithm_script, target_file, gorilla, topTF))

class HomeView(View):
    ''' This view serves serves as the index page of the website. '''

    # GET request (return index.html)
    def get(self, request, *args, **kwargs):
        form = UploadFileForm()
        context = {'navbar': 'tf_identification/includes/navbar.html',
                   'index': True,
                   'form': form,
                   'footer': 'tf_identification/includes/footer.html'}
        return render(request, "tf_identification/index.html", context)

    # POST file (file upload)
    def post(self, request, *args, **kwargs):

        form = UploadFileForm(request.POST, request.FILES)
        # check if form is valid
        if form.is_valid():
            # generate a unique ID for the process
            process_ID = get_random_string(length=20) + str(round(time.time()))
            # process the file
            go = 0
            if 'gorilla' in request.POST: go = 1
	    redirect("tfidentificationresults/" + process_ID + '/')
            handle_uploaded_file(request.FILES['file'],
                                 filename = process_ID,
                                 gorilla = go,
                                 topTF = request.POST['topTF']) # see definition above
            # return redirect("tfidentificationresults/" + process_ID + '/')
        else:
            pass
        # TODO: implement error handling if the form is not valid



class TFIdentificationResultView(View):
    ''' This view serves as the waiting and the result page for the TF identification algorithm.'''

    def get(self, request, *args, **kwargs):
        # get process ID
        if request.build_absolute_uri()[-1] == '/':
            process_ID = request.build_absolute_uri()[-31:-1]
        else:
            process_ID = request.build_absolute_uri()[-30:]

        # check if the results are ready for that link
        target_html = 'tf_identification/results/' + process_ID + '.html'
        if os.path.exists(BASE_DIR + '/tf_identification/templates/' + target_html):

            # the file exists, now return it
            while os.stat(BASE_DIR + '/tf_identification/templates/' + target_html).st_size < 1:
                # if the file is still being created
                time.sleep(0.1)
            context = {'navbar': 'tf_identification/includes/navbar.html',
                       'index': True,
                       'footer': 'tf_identification/includes/footer.html'}
            return render(request, target_html, context)

        else:
            # the file does not yet exist – show the loading page
            randomfacts = {'The human genome consists of 3.0 × 10^9 base pairs (3000 MB).',
                           'The largest cell in the human body is the female egg and the smallest is the male sperm.',
                           'The brain is much more active at night than during the day.',
                           'Your eyes are always the same size from birth (but your nose and ears never stop growing).',
                           "We have about the same number of hairs on our bodies as a chimpanzee, it's just that our hairs are useless, so fine they are almost invisible.",
                           }
            context = {'navbar': 'tf_identification/includes/navbar.html',
                       'index': True,
                       'refresh': 10,
                       'url': request.build_absolute_uri(),
                       'random_fact': random.sample(randomfacts, 1)[0],
                       'footer': 'tf_identification/includes/footer.html'}
            return render(request, "tf_identification/TF_identification_loading_page.html", context)




