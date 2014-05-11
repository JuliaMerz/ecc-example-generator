from app import app
import forms
from flask import render_template, request
import ecc_web

@app.route('/')
@app.route('/index')
def index():
    form = forms.ECCForm()
    return render_template('index.html', form=form)

@app.route('/ecc/', methods=['POST'])
def ecc():
    form = forms.ECCForm(request.form)
    if form.validate():
        try:
            output = ecc_web.run_web(int(form.n.data), int(form.k.data), int(form.mod.data), form.show_inv.data, form.show_rref.data, form.show_lagrange.data, form.show_orig.data, form.show_bw.data)
        except Exception:
            output = "ERROR: Please contact Sebastian Merz and tell him exactly what you did to cause this"
        return output
    return render_template("errors.html", form=form)
