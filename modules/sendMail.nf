process SENDMAIL {
    tag "sendMail"

    input:
    path(to_html)

    script:
    mail = [
        to: 'agletdinov@cmd.su',
        subject: 'Catch up',
        body: 'Report',
        attach: to_html
    ]

    sendMail(mail)
}