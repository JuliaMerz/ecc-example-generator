from app import app
import forms
from flask import render_template, request
import ecc_web

@app.route('/')
@app.route('/index')
def index():
    form = forms.ECCForm()
    return render_template('index.html', form=form)

@app.route('/ecc/', methods=['GET', 'POST'])
def ecc():
    print request.form
    form = forms.ECCForm(request.form)
    if form.validate():
        print form.k.data
        return ecc_web.run_web(int(form.n.data), int(form.k.data), int(form.mod.data), form.show_inv.data, form.show_rref.data, form.show_lagrange.data, form.show_orig.data, form.show_bw.data)
    return 'ERROR'
