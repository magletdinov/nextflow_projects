process SENDMAIL {
    tag "sendMail"

    input:
    path(to_html)

    script:
    mail = [
        to: 'agletdinov@cmd.su',
        subject: 'Catch up',
        body: 'Report',
        attach: "multiqc_report.html"
    ]

    sendMail(mail)
}

process SENDMAIL_PY {
    conda "/export/home/agletdinov/mambaforge/envs/reat"
    tag "sendMail"
    input:
    path(to_html)

    script:
    """
    #!/usr/bin/env python3
    from envelope import Envelope
    Envelope()\
        .from_("hpcreports@cmd.su")\
        .subject("[hpc_tasks] report")\
        .to("agletdinov@cmd.su")\
        .message("report")\
        .attach(name="multiqc_reporta.html", path="multiqc_report.html")\
        .smtp(host="email.cmd.su", port="25", user="hpcreports", password="tnXmk5k7KLy8", security="None", timeout=30, attempts=3, delay=3)\
        .send()
    """
}
