from allauth.account.models import EmailAddress, EmailConfirmation
from allauth.account.views import signup
from django.contrib.auth import get_user_model, login
from django.contrib.auth.models import AnonymousUser
from django.contrib.messages.middleware import MessageMiddleware
from django.contrib.sessions.middleware import SessionMiddleware
from django.test import TestCase, Client
from django.test.client import RequestFactory
from django.urls import reverse
from logging import getLogger

logger = getLogger(__name__)

User = get_user_model()
signup_url = reverse('account_signup')
login_url = reverse('account_login')
front_url = reverse('spectrum:index')

logger.debug(User)


class UserTestCase(TestCase):
    # signupページでユーザーを作成し、Eメール確認も済にする
    @classmethod
    def setUp(cls):
        cls.email = 'hoge@test.mail'
        cls.username = 'hoge-hoge'
        cls.password = 'jenh87dkIJFfwaosu42'
        # cls.res = cls.client.post(signup_url,
        #                             {'email': cls.email, 'username': cls.username,
        #                              'password1': cls.password, 'password2': cls.password})
        #
        # logger.debug(cls.res)
        # logger.debug(cls.res.url)

        cls.user_obj = User.objects.create(email=cls.email, username=cls.username, password=cls.password)
        logger.debug(cls.user_obj)

        cls.email_obj = EmailAddress.objects.create(user=cls.user_obj, email=cls.email,
                                                     primary=True, verified=True)

    # Try to register an existing account.
    def test_signup_existing_account(self):
        response = self.client.post(signup_url,
                                    {'email': self.email, 'username': self.username,
                                     'password1': self.password, 'password2': self.password})
        logger.debug(f'URL: {response.request["PATH_INFO"]}')
        logger.debug(f'Status code: {response.status_code}')

        self.assertEqual(response.request['PATH_INFO'], signup_url)
        self.assertEqual(response.status_code, 200)

        with open('./accounts/tests/sample_pages/test_signup_existing_account.html', 'w') as f:
            f.write(response.content.decode('utf-8'))

    def test_signup(self):
        response = self.client.post(signup_url,
                                    {'email': 'kanako.hirata@oist.jp', 'username': 'hoge_another',
                                     'password1': 'cQn5yt832Kfqlj1s', 'password2': 'cQn5yt832Kfqlj1s'})

        logger.debug(f'URL: {response.url}')
        logger.debug(f'Status code: {response.status_code}')

        self.assertEqual(response.url, '/accounts/confirm-email/')
        self.assertEqual(response.status_code, 302)  # Temporary redirect

        with open('./accounts/tests/sample_pages/test_signup.html', 'w') as f:
            f.write(response.content.decode('utf-8'))

    # ユーザーがきちんと作らているか？
    def test_single_user(self):
        self.assertEqual(User.objects.count(), 1)
        self.assertEqual(self.user_obj.email, self.email)

    # メールアドレスがverifiedになっているか？
    def test_emailaddress_verified(self):
        self.assertEqual(self.email_obj.verified, True)

    # ログインページでのpostが正しく動作するか
    def test_login_page(self):
        data = {"email": self.email, "password": self.password}
        response = self.client.post(login_url, data)
        self.assertEqual(response.status_code, 200)

    # ログアウトの際にログインページへリダイレクトされているか？
    def test_logout(self):
        res = self.client.get("/accounts/logout/")
        self.assertEqual(res.status_code, 302)
        self.assertEqual(res.url, front_url)
    #
    # # verifiedされているユーザーの場合はきちんとログインできているか？
    # def test_verified_user_login(self):
    #     c = Client()
    #     res_bool = c.login(email=self.email, password=self.password)
    #     self.assertEqual(res_bool, True)
    #
    # # verifiedされていないユーザーの場合はログインできない？
    # def test_not_verified_user_login(self):
    #     res = self.client.post(signup_url, {"email": "not-verifed@itc.tokyo", "password1": "somepassword",
    #                                         "password2": "somepassword"})
    #     self.assertEqual(res.status_code, 302)
    #     c = Client()
    #     res_bool = c.login(email="not-verified@itc.tokyo", password="somepassword")
    #     self.assertEqual(res_bool, False)